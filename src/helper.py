#!/usr/bin/python


import numpy as np

def sites_map(frag, xyz, sites):
    '''
    @description: map the binding sites to the 3D fragment

    @param:
            frag:   start, end; sorted
            xyz :   x, y, z
            sites:  position
    @return:
            sites_coord:    numpy array, coord of each site
            sites_id:       list, id of each site
    '''
    frag_idx = 0
    site_idx = 0
    sites_coord = []
    sites_id = []

    # sites center, sort
    sites_center = (sites[:,1] + sites[:,0]) // 2
    sites_center = np.sort(sites_center)

    while (frag_idx < len(frag)) and (site_idx < len(sites)):
        if sites_center[site_idx] <= frag[frag_idx,0]:
            site_idx += 1
        elif sites_center[site_idx] <= frag[frag_idx,1]:
            sites_coord.append(xyz[frag_idx])
            sites_id.append(frag_idx)
            site_idx += 1
        else:
            frag_idx += 1

    return np.array(sites_coord), sites_id


def degree_map(frag, region):
    '''
    @description: map the regions to 1D, get 1D fragment degrees

    @param:
            frag:   start, end; sorted
            region: position
    @return:
            frag_degree: numpy array, degree of each fragment
    '''
    frag_idx    = 0
    region_idx  = 0
    frag_degree = np.zeros(len(frag))

    # region center, sort
    region_center = (region[:,1] + region[:,0]) // 2
    region_center = np.sort(region_center)

    while (frag_idx < len(frag)) and (region_idx < len(region_center)):
        if region_center[region_idx] <= frag[frag_idx,0]:
            region_idx += 1
        elif region_center[region_idx] <= frag[frag_idx,1]:
            frag_degree[frag_idx] += 1
            region_idx += 1
        else:
            frag_idx += 1

    return frag_degree
