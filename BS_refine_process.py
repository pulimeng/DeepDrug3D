#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Visualize the effect of the refinement process:
1. Remove grid points close to the protein
2. Remove grid points outside the convex hull of the protein
3. Remove grid points that are not connected to the main pocket (top 3 largest clusters)

@author: limeng
"""

import numpy as np
import pandas as pd
from biopandas.mol2 import PandasMol2

from scipy.spatial.distance import cdist
import scipy
from scipy.spatial import Delaunay
from sklearn.cluster import DBSCAN

import string

# Parse the mol2 file
def parseMol2(pmol):
    atType = pmol.df['atom_type']
    coord = pmol.df[['x', 'y', 'z']]
    subset_name = pmol.df['subst_name']
    frames = pd.concat([atType,coord,subset_name], axis = 1)
    return frames
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
    tree = scipy.spatial.cKDTree(data)
    mask = tree.query_ball_point(center,1.01*r)
    points_in_sphere = data[mask]
    return points_in_sphere
# make line for pdb file
def makeLine(coord, potE, k, num_chr_dict):
    atom = 'ATOM'
    atom_sn = str(k+1)
    atom_name = 'D1'
    hundred, one = np.divmod(k,676)
    ten, one = np.divmod(k-hundred*676,26)
    res_name = num_chr_dict[hundred] + num_chr_dict[ten] + num_chr_dict[one]
    E = potE
    x = '{:.3f}'.format(round(coord[0], 3))
    y = '{:.3f}'.format(round(coord[1], 3))
    z = '{:.3f}'.format(round(coord[2], 3))
    EE = '{:.2f}'.format(round(E,2))
    string = atom + ' '*2 + '{:>5}'.format(atom_sn) + ' ' + '{:4}'.format(atom_name) + ' ' \
            + '{:>3}'.format(res_name) + ' '*2 + '   1' + ' '*4 + '{:>8}'.format(x) + '{:>8}'.format(y) + '{:>8}'.format(z) \
            + '{:>6}'.format('1.00') + '{:>6}'.format(EE) + ' '*8 + '\n'
    return string
# test if a point is inside a convex hull
def in_hull(p, hull):
    if not isinstance(hull,Delaunay):
        hull = Delaunay(hull)
    return hull.find_simplex(p)>=0
# load the data
protein_name = '1a2sA'
protein_path = '/home/limeng/Desktop/pocket_similarity/binding_data/all_protein-mol2/' + protein_name + '.mol2'
ligand_name = '1a2sA00'
ligand_path = '/home/limeng/Desktop/pocket_similarity/binding_data/all_ligand-mol2/' + ligand_name + '.mol2'
protein = PandasMol2().read_mol2(protein_path)
ligand = PandasMol2().read_mol2(ligand_path)
pt = parseMol2(protein)
ld = parseMol2(ligand)
pt_coords = np.array([pt['x'],pt['y'],pt['z']]).T
ld_coords = np.array([ld['x'],ld['y'],ld['z']]).T
ld_center = np.mean(ld_coords, axis = 0)
num_range = np.linspace(0,25,26, dtype = int)
chr_range = list(string.ascii_uppercase[:27])
num_chr_dict = dict(zip(num_range,chr_range))
# original sphere BS
BS_original = sGrid(ld_center, 15, 31)
BS_original = np.concatenate((BS_original,np.zeros((len(BS_original),1))), axis = 1)
with open (protein_name + '_site_original.pdb','w') as in_strm:
    for k in xrange(len(BS_original)):
        temp_coords = BS_original[k,:]
        cl = temp_coords[-1]
        temp_string = makeLine(temp_coords, cl, k, num_chr_dict)
        in_strm.write(temp_string)
print 'Number of points is ' + str(len(BS_original))
# Remove points too close to the protein (less than 2 A)
dist = cdist(BS_original[:,0:3], pt_coords, 'euclidean')
BS_2close1 = []
for i in xrange(len(dist)):
    temp = BS_original[i,:]
    if np.any(dist[i,:] < 2.1):
        temp[-1] = 1
        BS_2close1.append(temp)
    else:
        temp[-1] = 0
        BS_2close1.append(temp)
BS_2close1 = np.array(BS_2close1)
with open (protein_name + '_site_2close1.pdb','w') as in_strm:
    for k in xrange(len(BS_2close1)):
        temp_coords = BS_2close1[k,:]
        cl = temp_coords[-1]
        temp_string = makeLine(temp_coords, cl, k, num_chr_dict)
        in_strm.write(temp_string)
BS_2close2 = []
for i in xrange(len(dist)):
    temp = BS_original[i,:]
    if np.any(dist[i,:] < 2.1):
        continue
    else:
        BS_2close2.append(temp)
BS_2close2 = np.array(BS_2close2)
with open (protein_name + '_site_2close2.pdb','w') as in_strm:
    for k in xrange(len(BS_2close2)):
        temp_coords = BS_2close2[k,:]
        cl = temp_coords[-1]
        temp_string = makeLine(temp_coords, cl, k, num_chr_dict)
        in_strm.write(temp_string)
print 'Number of points after removing points close to the protein is ' + str(len(BS_2close1))
# Remove points outside the convex hull
hull_labels = in_hull(BS_2close2[:,0:3], pt_coords)
BS_inhull1 = BS_2close2[hull_labels]
with open (protein_name + '_site_inhull1.pdb','w') as in_strm:
    for k in xrange(len(BS_inhull1)):
        temp = BS_inhull1[k,:]
        cl = temp[3]
        temp_string = makeLine(temp, cl, k, num_chr_dict)
        in_strm.write(temp_string)
BS_inhull2 = []
for k in xrange(len(hull_labels)):
    temp = BS_2close2[k,:]
    if hull_labels[k]:
        temp[-1] = 1
        BS_inhull2.append(temp)
    else:
        temp[-1] = 0
        BS_inhull2.append(temp)
BS_inhull2 = np.array(BS_inhull2)
with open (protein_name + '_site_inhull2.pdb','w') as in_strm:
    for k in xrange(len(BS_inhull2)):
        temp = BS_inhull2[k,:]
        cl = temp[3]
        temp_string = makeLine(temp, cl, k, num_chr_dict)
        in_strm.write(temp_string)
print 'Number of points after removing points outside the convex hull of the protein is ' + str(len(BS_inhull1))
# remove isolated grid points
iso_dist = cdist(BS_inhull1[:,0:3],BS_inhull1[:,0:3])
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
iso_site = []
cnt = 0
for j in xrange(len(final_label)):
    if final_label[j] == 1:
        temp = np.append(BS_inhull1[j,0:3],np.array(1))
        iso_site.append(temp)
        cnt += 1
    else:
        temp = np.append(BS_inhull1[j,0:3],np.array(0))
        iso_site.append(temp)
iso_site = np.array(iso_site)
with open (protein_name + '_site_cc2.pdb','w') as in_strm:
    for k in xrange(len(iso_site)):
        temp_c = iso_site[k,:]
        potE_c = temp_c[3]
        temp_string = makeLine(temp_c, potE_c, k, num_chr_dict)
        in_strm.write(temp_string)
print 'Number of points after removing points not connected to the main site is ' + str(cnt)
#iso_site = data_c
#with open (protein_name + '_site_centered.pdb','w') as in_strm:
#    for k in xrange(len(iso_site)):
#        temp_c = iso_site[k,:]
#        potE_c = temp_c[3]
#        temp_string = makeLine(temp_c, potE_c, k, num_chr_dict)
#        in_strm.write(temp_string)
#print 'Number of points after removing points not connected to the main site is ' + str(cnt)
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
    print len(inside_site)
    # remove any grid points outside the convex hull
    in_bool = in_hull(inside_site[:,0:3], protein_coords)
    hull_site = inside_site[in_bool]
    print len(hull_site)
    # remove isolated grid points
    iso_dist = cdist(hull_site[:,0:3],hull_site[:,0:3])
    labels = DBSCAN(eps = 1.414, min_samples = 3, metric = 'precomputed').fit_predict(iso_dist)
    unique, count = np.unique(labels, return_counts = True)
    sorted_label = [x for _,x in sorted(zip(count,unique))]
    sorted_label = np.array(sorted_label)
    null_index = np.argwhere(sorted_label == -1)
    cluster_labels = np.delete(sorted_label, null_index)
    save_labels = np.flip(cluster_labels, axis = 0)[0] #TODO change this to the single largest CC
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

aa = site_refine(BS_original,pt_coords)