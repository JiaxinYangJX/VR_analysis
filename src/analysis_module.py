#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Author: Jiaxin Yang
Date: 2021-04-22 13:52:51
LastEditTime: 2021-04-22 18:06:32
LastEditors: Jiaxin Yang
Description: 
FilePath: \undefinedd:\Program Files (x86)\vs_python\VR_analysis\src\analysis_module.py
'''

import pandas as pd
import numpy as np



def spatial_hub(sturcture, sites):
    '''
    @description: generate 3D spatial hubs of specific sites
    @param:
            structure:  3D structurs
            sites:      a specific binding sites, protein, DNase, gv
    @return:
    '''
    # cumulative counts
    # counts difference matrix
    # density matrix
    # find out

    return output


def interaction_hub(sturcture, links):
    '''
    @description: generate hubs with high degree of interaction
    @param:
            structure:  3D structurs
            links:      interaction data, eQTL
    @return:
    '''
    
    return output


def loop_3d(structure, size, threshold):
    '''
    @description: get the chromatin loops 
    @param:
            structure:  3D structurs
            size:       size of loops
            threshold:  distance threshold
    @return:
    '''
    # calcualte pairwise distance
    # find out the points with small 3D distance
    # filter the size

    return output


def neighbor_3D(structure, distance):
    '''
    @description: get the closes module
    @param:
            structure:  3D structurs
            distance:   distance threshold
    @return:
    '''
    # seems similar with loop
    return output