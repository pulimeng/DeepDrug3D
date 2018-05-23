#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This script is for the training of our model to classify binidng site ATP- or Heme-
binding.
The input should be the voxel representation of the pocket of the protein and the
labels.

@author: Limeng Pu
"""

import sys
import os
import argparse

import numpy as np

from deepdrug3d import DeepDrug3DBuilder

from keras import callbacks
from keras.optimizers import Adam
from keras.utils import np_utils

def argdet():
    if len(sys.argv) < 13:
        print('Check number of input arguement!')
        exit()
    elif len(sys.argv) == 13:
        print('Save model as deepdrug3d.h5')
        args = myargs()
        return args
    elif len(sys.argv) == 15:
        print('Save model to provided location')
        args = myargs()
        return args
    else:
        print('Cannot recognize the inputs!')
        exit()

def myargs():
    parser = argparse.ArgumentParser()                                              
    parser.add_argument('--alist', required = True, help = 
                        'location of list contains names of proteins binds with ATP')
    parser.add_argument('--hlist', required = False, help = 
                        'location of list contains names of proteins binds with Heme')
    parser.add_argument('--vfolder', required = True, help = 'folder for the voxel data')
    parser.add_argument('--bs', required = True, help = 'batch size')
    parser.add_argument('--lr', required = True, help = 'initial learning rate')
    parser.add_argument('--epoch', required = True, help = 'number of epochs for taining')
    parser.add_argument('--output', required = False, help = 'location for the model to be saved')
    args = parser.parse_args()
    return args
    
def train_deepdrug(atp_list, heme_list, voxel_folder, batch_size, lr, epoch, output):
    mdl = DeepDrug3DBuilder.build()
    adam = Adam(lr=lr, beta_1=0.9, beta_2=0.999, epsilon=None, decay=0.0, amsgrad=False)
    
    # We add metrics to get more results you want to see
    mdl.compile(optimizer=adam, loss='binary_crossentropy', metrics=['accuracy'])
    
    # load the data
    atps = []
    with open(atp_list) as atp_in:
        for line in atp_in.readlines():
            temp = line.replace(' ','').replace('\n','')
            atps.append(temp)
    hemes = []
    with open(heme_list) as heme_in:
        for line in heme_in.readlines():
            temp = line.replace(' ','').replace('\n','')
            hemes.append(temp)
    # conver data into a single matrix
    atp_len = len(atps)
    heme_len = len(hemes)
    L = atp_len + heme_len
    voxel = np.zeros(shape = (L, 14, 32, 32, 32),
        dtype = np.float64)
    label = np.zeros(shape = (L,), dtype = int)
    cnt = 0
    print('...Loading the data')
    for filename in os.listdir(voxel_folder):
        protein_name = filename[0:-4]
        full_path = voxel_folder + '/' + filename
        temp = np.load(full_path)
        voxel[cnt,:] = temp
        if protein_name in atps:
                label[cnt] = 0
        elif protein_name in hemes:
                label[cnt] = 1
        else:
            print protein_name
            print 'Something is wrong...'
            break
        cnt += 1
        
    y = np_utils.to_categorical(label, num_classes = 2)
    # callback function for model checking
    tfCallBack = callbacks.TensorBoard(log_dir='./graph', histogram_freq = 0, batch_size=batch_size, write_graph=True, 
    					write_grads=False, write_images=True, embeddings_freq=0, embeddings_layer_names=None, embeddings_metadata=None)
    mdl.fit(voxel, y, epochs = epoch, batch_size = batch_size, 
              shuffle = True, callbacks = [tfCallBack], verbose = 2)
    # save the model
    if output == None:
        mdl.save('deepdrug3d.h5')
    else:
        mdl.save(output)
    
if __name__ == "__main__":
    args = argdet()
    train_deepdrug(args.alist, args.hlist, args.vfolder, args.batch_size, args.lr, args.epoch, args.output)
