#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Generate the binding grid and calculate the DFIRE potentials

@author: Limeng Pu
"""

import numpy as np
import scipy.spatial as sp
import subprocess
import os

import random
import string

import time

import pandas as pd
from scipy.spatial.distance import cdist
from scipy.spatial import Delaunay
from sklearn.cluster import DBSCAN

# Generate sphere grid points
def sGrid(center, r, N):
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
# test if a point is inside a convex hull
def in_hull(p, hull):
    if not isinstance(hull,Delaunay):
        hull = Delaunay(hull)
    return hull.find_simplex(p)>=0
# binding site refinement
def site_refine(site,protein_coords):
    # distance matrix for the removal of the grid points that are too close (<= 2 A) to any protein atoms
    dist = cdist(site[:,0:3], protein_coords, 'euclidean')
    inside_site = []
    for i in xrange(len(dist)):
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
    for k in xrange(len(labels)):
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
The following calculates the DFIRE potentials using the dligand program proivded in the DFIRE paper
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
        a = out[0].replace('\n','')
        b = float(a)
        append(b)
    Es = np.array(Es)
    os.remove(temp_filename)
    return Es
def minmax_scale(X,axis = 0):
    X_std = (X - X.min(axis=0)) / (X.max(axis=0) - X.min(axis=0))
    X_scaled = X_std * (X.max - X.min) + X.min
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
    def build(protein_coords, protein_path, r, N):
        """
        Input: protein coordinates, path to the pdb file of the protein, radius, number of points along the radius.
        Output: dataframe of the binding grid, including coordinates and potentials for different atom types.
        """
        print 'The radius of the binding grid is: ' + str(r)
        print 'The number of points along the diameter is: ' + str(N)
        binding_site = sGrid(np.array([0,0,0]),r,N)
        new_site = site_refine(binding_site, protein_coords)
        print 'The number of points in the refined binding set is ' + str(len(new_site))
        ss = time.time()
        print '... Computation of the binding site potential energy started ...'
        total_potE = potEnergy(new_site, 'dummy_mol2.mol2', protein_path)
        print 'The total time of binding site potential energy computation is: ' + str(time.time() - ss) + ' seconds'
        df1 = pd.DataFrame(total_potE['loc'], columns = ['x','y','z'])
        df2 = pd.DataFrame(total_potE['potE'], columns = ['C.2','C.3','C.ar','F','N.am','N.2','O.co2','N.ar','S.3','O.2','O.3','N.4','P.3','N.pl3'])
        frames = [df1,df2]
        df = pd.concat(frames, axis = 1)
        return  df