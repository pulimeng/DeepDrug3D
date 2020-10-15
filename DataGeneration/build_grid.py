import argparse

import numpy as np
import scipy.spatial as sp
import subprocess
import os
import os.path as osp

import random
import string

import time

import pandas as pd
from openbabel import pybel
from biopandas.pdb import PandasPdb
from scipy.spatial.distance import cdist
from scipy.spatial import Delaunay
from sklearn.cluster import DBSCAN

"""
The following functions process the input pdb files by moving the center of pocket to [0,0,0]
and align the protein to the principle axes of the pocket
"""
def normalize(v):
    """
    vector normalization
    """
    norm = np.linalg.norm(v)
    if norm == 0: 
       return v
    return v / norm

def vrrotvec(a,b):
    """
    Function to rotate one vector to another, inspired by
    vrrotvec.m in MATLAB
    """
    a = normalize(a)
    b = normalize(b)
    ax = normalize(np.cross(a,b))
    angle = np.arccos(np.minimum(np.dot(a,b),[1]))
    if not np.any(ax):
        absa = np.abs(a)
        mind = np.argmin(absa)
        c = np.zeros((1,3))
        c[mind] = 0
        ax = normalize(np.cross(a,c))
    r = np.concatenate((ax,angle))
    return r

def vrrotvec2mat(r):
    """ 
    Convert the axis-angle representation to the matrix representation of the 
    rotation 
    """
    s = np.sin(r[3])
    c = np.cos(r[3])
    t = 1 - c
    
    n = normalize(r[0:3])
    
    x = n[0]
    y = n[1]
    z = n[2]
    
    m = np.array(
     [[t*x*x + c,    t*x*y - s*z,  t*x*z + s*y],
     [t*x*y + s*z,  t*y*y + c,    t*y*z - s*x],
     [t*x*z - s*y,  t*y*z + s*x,  t*z*z + c]]
    )
    return m

def coords_transform(protein_coords, pocket_center, pocket_coords):
    """
    Transform the protein coordinates so that the pocket is centered at [0,0,0]
    and align the protein coordinates according to the principle axes of the pocket
    """
    pocket_coords = pocket_coords - pocket_center # center the pocket to 0,0,0
    protein_coords = protein_coords - pocket_center # center the protein according to the pocket center
    
    inertia = np.cov(pocket_coords.T)
    e_values, e_vectors = np.linalg.eig(inertia)
    sorted_index = np.argsort(e_values)[::-1]
    sorted_vectors = e_vectors[:,sorted_index]
    # Align the first principal axes to the X-axes
    rx = vrrotvec(np.array([1,0,0]),sorted_vectors[:,0])
    mx = vrrotvec2mat(rx)
    pa1 = np.matmul(mx.T,sorted_vectors)
    # Align the second principal axes to the Y-axes
    ry = vrrotvec(np.array([0,1,0]),pa1[:,1])
    my = vrrotvec2mat(ry)
    transformation_matrix = np.matmul(my.T,mx.T)
    # transform the protein coordinates to the center of the pocket and align with the principal
    # axes with the pocket
    transformed_coords = (np.matmul(transformation_matrix,protein_coords.T)).T
    return transformed_coords

"""
The following functions generate and refine the binding pocket grid
"""
def sGrid(center, r, N):
    """
    Generate spherical grid points at the center provided
    """
    center = np.array(center)
    x = np.linspace(center[0]-r,center[0]+r,N)
    y = np.linspace(center[1]-r,center[1]+r,N)
    z = np.linspace(center[2]-r,center[2]+r,N)
    #Generate grid of points
    X,Y,Z = np.meshgrid(x,y,z)
    data = np.vstack((X.ravel(),Y.ravel(),Z.ravel())).T
    # indexing the interior points
    tree = sp.cKDTree(data)
    mask = tree.query_ball_point(center,1.01*r)
    points_in_sphere = data[mask]
    return points_in_sphere

def in_hull(p, hull):
    """
    Test if a point is inside a convex hull
    """
    if not isinstance(hull,Delaunay):
        hull = Delaunay(hull)
    return hull.find_simplex(p)>=0

def site_refine(site, protein_coords):
    """
    Binding site refinement
    """
    # distance matrix for the removal of the grid points that are too close (<= 2 A) to any protein atoms
    dist = cdist(site[:,0:3], protein_coords, 'euclidean')
    inside_site = []
    for i in range(len(dist)):
        if np.any(dist[i,:] < 2.1):
            continue
        else:
            inside_site.append(site[i,:])
    inside_site = np.array(inside_site)
    # remove any grid points outside the convex hull
    in_bool = in_hull(inside_site[:,0:3], protein_coords)
    hull_site = inside_site[in_bool]
    # remove isolated grid points
    iso_dist = cdist(hull_site[:,0:3],hull_site[:,0:3])
    labels = DBSCAN(eps = 1.414, min_samples = 3, metric = 'precomputed').fit_predict(iso_dist)
    unique, count = np.unique(labels, return_counts = True)
    sorted_label = [x for _,x in sorted(zip(count,unique))]
    sorted_label = np.array(sorted_label)
    null_index = np.argwhere(sorted_label == -1)
    cluster_labels = np.delete(sorted_label, null_index)
    save_labels = np.flip(cluster_labels, axis = 0)[0]
    final_label = np.zeros(labels.shape)
    for k in range(len(labels)):
        if labels[k] == save_labels:
            final_label[k] = 1
        else:
            continue
    final_label = np.array(final_label, dtype = bool)
    # potential energy normalization
    iso_site = hull_site[final_label]
    return iso_site

""" 
The following functions create a new dummy mol2 file for the DFIRE calculation 
"""
# replace the coordinates in the original string with new coordinates
def replace_coord(original_string, new_coord):
    temp = '{:>8}  {:>8}  {:>8}'.format(new_coord[0],new_coord[1],new_coord[2])
    new_string = original_string.replace(' 50.0000   51.0000   52.0000',temp)
    return new_string

# replace the atom type in the original string with the new atom type
def replace_type(original_string, new_type):
    temp = '{:6}'.format(new_type)
    new_string = original_string.replace('N.3   ',temp)
    return new_string

# replace the residue type with new residue type
def replace_res(original_string, new_res):
    temp = '{:6}'.format(new_res)
    new_string = original_string.replace('VAL1  ',temp)
    return new_string

"""
The following functions calculate the DFIRE potentials using the dligand program proivded in the DFIRE paper
"""
# UISNG THE DFIRE FUNCTION
def single_potEnergy(loc1, ld_type_list, mol2_in_string, protein_file):
    temp_loc = loc1.round(4)
    Es = []
    append = Es.append
    r1 = replace_coord(mol2_in_string, temp_loc)
    random_string = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(11))
    temp_filename = '/dev/shm/' + random_string +'.mol2' # TODO: this controls the place to generate the temporary mol2 file
    for item in ld_type_list:
        rrr = replace_type(r1, item)
        f = open(temp_filename,'w')
        f.write(rrr)
        f.close()
        child = subprocess.Popen(['./dligand-linux', temp_filename, protein_file],stdout=subprocess.PIPE)
        child.wait()
        out = child.communicate()
        out = out[0].decode("utf-8") 
        a = out.replace('\n','')
        b = float(a)
        append(b)
    Es = np.array(Es)
    os.remove(temp_filename)
    return Es

def minmax_scale(X,axis = 0):
    X_std = (X - X.min(axis=0)) / (X.max(axis=0) - X.min(axis=0))
    X_scaled = X_std * (X.max() - X.min()) + X.min()
    return X_scaled

def potEnergy(binding_site, mol2_in_path, protein_file):
    ld_type_list = ['C.2','C.3','C.ar','F','N.am','N.2','O.co2','N.ar','S.3','O.2','O.3','N.4','P.3','N.pl3']
    total_potE = {'loc':[],'potE':[]}
    mol2_in_file = open(mol2_in_path)
    mol2_in_string = mol2_in_file.read()
    potEs = np.array([single_potEnergy(loc1, ld_type_list, mol2_in_string, protein_file) for loc1 in binding_site])
    total_potE['potE'] = minmax_scale(potEs, axis = 0)
    total_potE['loc'] = binding_site
    return total_potE

# main function
class Grid3DBuilder(object):
    """ Given an align protein, generate the binding grid 
    and calculate the DFIRE potentials """
    @staticmethod
    def build(pdb_path, aux_input_path, r, N, output_folder, shape):
        """
        Input: protein coordinates, path to the pdb file of the protein, radius, number of points along the radius.
        Output: dataframe of the binding grid, including coordinates and potentials for different atom types.
        """
        # Parse the pdb file and the auxilary input file
        ppdb = PandasPdb().read_pdb(pdb_path)
        protein_df = ppdb.df['ATOM']
        content = []
        with open(aux_input_path) as in_strm:
            for line in in_strm.readlines():
                l = line.replace('\n','')
                idx = l.index(':')
                content.append(l[idx+1:None])
        resi =  [int(x) for x in content[0].split(' ')]
        if len(content[1]) != 0:
            pocket_center = np.array([int(x) for x in content[1].split(' ')])
            print('Center provided as {:.2f} {:.2f} {:.2f}'.format(pocket_center[0], pocket_center[1], pocket_center[2]))
        else:
            print('No center is provided')
            pocket_df = protein_df[protein_df['residue_number'].isin(resi)]
            pocket_coords = np.array([pocket_df['x_coord'], pocket_df['y_coord'], pocket_df['z_coord']]).T
            pocket_center = np.mean(pocket_coords, axis = 0)
            print('Center calculated as {:.2f} {:.2f} {:.2f}'.format(pocket_center[0], pocket_center[1], pocket_center[2]))
        protein_coords = np.array([protein_df['x_coord'], protein_df['y_coord'], protein_df['z_coord']]).T
        transformed_coords = coords_transform(protein_coords, pocket_center, pocket_coords)
        # Generate a new pdb file with transformed coordinates
        ppdb.df['ATOM']['x_coord'] = transformed_coords[:,0]
        ppdb.df['ATOM']['y_coord'] = transformed_coords[:,1]
        ppdb.df['ATOM']['z_coord'] = transformed_coords[:,2]
        output_trans_pdb_path = osp.join(output_folder, pdb_path[:-4] + '_transformed.pdb')
        print('Saving the binding pocket aligned pdb file to: {}.pdb'.format(pdb_path[:-4]))
        ppdb.to_pdb(output_trans_pdb_path)
        output_trans_mol2_path = osp.join(output_folder, pdb_path[:-4] + '_transformed.mol2')
        print('Saving the binding pocket aligned mol2 file to: {}.mol2'.format(pdb_path[:-4]))
        mol = next(pybel.readfile('pdb',output_trans_pdb_path))
        mol.write('mol2',output_trans_mol2_path, overwrite=True)
        
        # Grid generation and DFIRE potential calculation
        print('The radius of the binding grid is: {}'.format(r))
        print('The number of points along the diameter is: {}'.format(N))
        binding_site = sGrid(np.array([0,0,0]), r, N)
        new_site = site_refine(binding_site, transformed_coords)
        print('The number of points in the refined binding set is {}'.format(len(new_site)))
        ss = time.time()
        if shape == True:
            print('Output only shape of the binidng grid')
            total_potE = new_site
            df = pd.DataFrame(total_potE, columns = ['x','y','z'])
            df.to_csv(osp.join(output_folder, pdb_path[:-4] + '.grid'), index=False)
        else:
            print('Calculating of the binding site potential energy')
            total_potE = potEnergy(new_site, 'dummy_mol2.mol2', output_trans_mol2_path)
            print('The total time of binding site potential energy computation is: {:.4f} seconds'.format(time.time() - ss))
            df1 = pd.DataFrame(total_potE['loc'], columns = ['x','y','z'])
            df2 = pd.DataFrame(total_potE['potE'], columns = ['C.2','C.3','C.ar','F','N.am','N.2','O.co2','N.ar','S.3','O.2','O.3','N.4','P.3','N.pl3'])
            frames = [df1,df2]
            df = pd.concat(frames, axis=1)
            df.to_csv(osp.join(output_folder, pdb_path[:-4] + '.grid'), index=False)
        return  df
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--f', type=str, required=True, help='Input pdb file name')
    parser.add_argument('--a', type=str, required=True, help='Input auxilary file name')
    parser.add_argument('--o', type=str, required=True, help='Output folder name')
    parser.add_argument('--r', type=int, required=True, help='Grid radius')
    parser.add_argument('--n', type=int, required=True, help='Number of points along diameter')
    parser.add_argument('--s', dest='shape', action='store_true', help='Return shape of grid only')
    parser.add_argument('--p', dest='shape', action='store_false')
    parser.set_defaults(shape=False)

    opt = parser.parse_args()
    Grid3DBuilder().build(opt.f, opt.a, opt.r, opt.n, opt.o, opt.shape)