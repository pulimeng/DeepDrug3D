#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This script contains the model of the proposed DeepDrug3D network. It is written using
keras with tensorflow as backend.

@author: Limeng Pu
"""

import os
os.environ['KERAS_BACKEND']='tensorflow'

from keras.models import Sequential
from keras.layers import Dense, Activation, Convolution3D, MaxPooling3D, Dropout, Flatten
from keras.layers.advanced_activations import LeakyReLU
from keras import backend as K

# For reproductivity
#seed = 12306
#np.random.seed(seed)

class DeepDrug3DBuilder(object):
    """DeepDrug3D network
    Ouput: a model takes a 5D tensor as input and outpus probabilities of the pocket binds to 
    ATP and Heme.
    """
    @staticmethod
    def build():
        model = Sequential()
        # Conv layer 1
        model.add(Convolution3D(
            input_shape = (14,32,32,32),
            filters=64,
            kernel_size=5,
            padding='valid',     # Padding method
            data_format='channels_first',
        ))
        model.add(LeakyReLU(alpha = 0.1))
        # Dropout 1
        model.add(Dropout(0.2))
        # Conv layer 2
        model.add(Convolution3D(
            filters=64,
            kernel_size=3,
            padding='valid',     # Padding method
            data_format='channels_first',
        ))
        model.add(LeakyReLU(alpha = 0.1))
        # Maxpooling 1
        model.add(MaxPooling3D(
            pool_size=(2,2,2),
            strides=None,
            padding='valid',    # Padding method
            data_format='channels_first'
        ))
        # Dropout 2
        model.add(Dropout(0.4))
        # FC 1
        model.add(Flatten())
        model.add(Dense(128)) # TODO changed to 64 for the CAM
        model.add(LeakyReLU(alpha = 0.1))
        # Dropout 3
        model.add(Dropout(0.4))
        # Fully connected layer 2 to shape (2) for 2 classes
        model.add(Dense(2))
        model.add(Activation('softmax'))
        return model