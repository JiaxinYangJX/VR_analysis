#!/usr/bin/python

import pandas as pd
import numpy as np
from helper import *
from sklearn.cluster import KMeans, AgglomerativeClustering
from collections import Counter



def read_structure(xyz_path, frag_id_path):
    '''
    @description: read structure

    @param:
            xyz_path:       path to xyz coord
            frag_id_path:   path to fragment 1D       
    @return:
            frag:   numpy array, fragment 1D position
            xyz:    numpy array, xyz
    '''    
    xyz     = pd.read_csv(xyz_path, header=None, sep='\t')
    frag_id = pd.read_csv(frag_id_path, header=None, sep='\t')
    
    xyz   = xyz.iloc[:,1:].values
    start = frag_id.iloc[:,0].values
    start = np.expand_dims(start, axis=1) 
    end   = start + 5000
    frag  = np.concatenate([start, end], axis=1)

    return frag, xyz


def read_sites(sites_path, chr_id):
    '''
    @description: read binding sites dataset

    @param:
            sites_path: path to the binding sites data
            chr_id:     chr id
    @return:
            sites: numpy array, binding sites
    '''
    sites = pd.read_csv(sites_path, header=None, sep='\t')
    sites = sites[sites.iloc[:,0]==chr_id].iloc[:,1:3].values

    return sites 


def read_links(links_path, chr_id):
    '''
    @description: read chromatin interaction dataset

    @param:
            links_path: path to the links: chr, frag1_start, end, frag2_start, end
            chr_id:     chr id
    @return:
            links:      numpy array, fragment links
    '''
    links = pd.read_csv(links_path, header=None, sep='\t')
    links = links[links.iloc[:,0]==chr_id].iloc[:,1:5].values
    return links


'''
def spatial_hub(frag, xyz, sites, cluster_size=40):
    ''
    @description: generate 3D spatial hubs of specific sites
    @param:
            structure:  3D structurs: start end x y z
            sites:      a specific binding sites, protein, DNase, gv
    @return:
    ''
    # map the sites into 3D structure
    sites_coord, sites_id = sites_map(frag, xyz, sites)
    
    # k-means clustering
    ncluster = len(sites_coord) // cluster_size
    my_kmeans = KMeans(n_clusters=ncluster).fit(sites_coord)
    
    # filter out the broad clusers
    for i in range(ncluster):
        dist = sites_coord[my_kmeans.labels_==i] - my_kmeans.cluster_centers_[i]
        print(np.sum(np.sum(dist**2, axis=1) < 10))



    return output
'''

def spatial_hub_hiera(frag, xyz, sites, dist_thres=3, cluster_size_thres=0.5):
    '''
    @description: generate 3D spatial hubs of specific sites

    @param:
            frag:       frag_id
            xyz:        xyz
            sites:      a specific binding sites, protein, DNase, gv
            dist_thres: distance threshold in hierachial clustering
            cluster_size_thres: only clusters with top 0.5 sizes
    @return:
            group_list: list, contains the frag_id in each 3D hub
    '''
    # 1. map the sites into 3D structure
    sites_coord, sites_id = sites_map(frag, xyz, sites)
    
    # 2. hierachical cluster
    my_hiera = AgglomerativeClustering(distance_threshold = dist_thres).fit(sites_coord)

    # 3. only keep the cluster with enough fragments, default: top 50%
    cluster_counter = Counter(my_hiera.labels_)
    size_thres = np.quantile(np.array(list(cluster_counter.items()))[:,1],q=cluster_size_thres)
    group_list = []
    for label, count in cluster_counter.most_common():
        if count > size_thres:
            group_list.append(sites_coord[my_hiera.labels_==label,])
 
    return group_list




def interaction_hub(frag, xyz, links, q_quantile = 0.99):
    '''
    @description: generate hubs with high degree of interaction

    @param:
            frag:       frag_id
            xyz:        xyz
            links:      frag-frag links
            q_quantile: top 0.99 dense degree
    @return:
            group_list: numpy array, contains the start and end id of hubs
    '''
    # 1. links to 2 regions
    region_1 = links[:,0:2]
    region_2 = links[:,2:4]
    region   = np.concatenate([region_1, region_2], axis=0)

    # 2. map to 1d, get degree
    frag_degree = degree_map(frag, region)

    # 3. cumulative increase
    cum_degree = np.cumsum(frag_degree)

    # 4. find the dense 1D region
    size = 5
    degree_list = []
    for p in range(frag.shape[0]-size):
        degree_list.append(cum_degree[p+size]-cum_degree[p])
    degree_list = np.array(degree_list)
    
    # find the high degree regions
    thres = np.quantile(degree_list, q = q_quantile)
    high_region_start = np.where(degree_list > thres)[0] # high range: (p,p+size]
    idx = 0
    start_idx = high_region_start[0] + 1 # [p+1,p+size]
    # merge the region
    group_list = []
    while idx < len(high_region_start)-1:
        if (high_region_start[idx] + size) >= high_region_start[idx+1]:
            # overlap
            idx += 1
        else: # save
            group_list.append([start_idx, high_region_start[idx]+size])
            start_idx = high_region_start[idx+1] + 1
            idx += 1
    group_list.append([start_idx, high_region_start[idx]+size]) # add last
    
    return np.array(group_list)



def loop_3d(frag, xyz, scale = 100000, resolution = 5000, q_quantile=0.002):
    '''
    @description: get the chromatin loops

    @param:
            frag:       frag_id
            xyz:        xyz
            scale:      loop scale
            resolution: resolution of the structure
            q_quantile: top 0.002 closest 
    @return:
            loop_list:  numpy array, contains the start and end id of loops
    '''
    # 1. find the 1) distance between two fragment
    size = scale // resolution
    dist_list = []
    for p in range(frag.shape[0]-size+1):
        dist_tmp = np.linalg.norm(xyz[p] - xyz[p+size-1])
        dist_list.append(dist_tmp)
    dist_list = np.array(dist_list)

    # 2. find the loop 
    thres = np.quantile(dist_list, q = q_quantile)
    close_loop_start = np.where(dist_list < thres)[0] # range: [p,p+size]
    
    # 3. merge the loop
    idx = 0
    start_idx = close_loop_start[0] 
    loop_list = []
    while idx < len(close_loop_start)-1:
        if (close_loop_start[idx] + size) >= close_loop_start[idx+1]:
            # overlap
            idx += 1
        else: # save
            loop_list.append([start_idx, start_idx + size])
            start_idx = close_loop_start[idx+1]
            idx += 1
    loop_list.append([start_idx, start_idx+size]) # add last
    
    return np.array(loop_list)


def main():
    xyz_path     = '../data/structure/chr1_1502144569709.xyz.txt'
    frag_id_path = '../data/structure/chr1_coordinate_mapping.txt'
    sites_path   = '../data/binding/ENCSR000EMT_rep2_1_se_bwa_biorep_filtered_peaks.bed'
    links_path   = '../data/links/GM_link.txt'
    chr_id       = 'chr1'

    frag, xyz = read_structure(xyz_path, frag_id_path)
    sites = read_sites(sites_path, chr_id)
    links = read_links(links_path, chr_id)
    
    hub_3d = spatial_hub_hiera(frag, xyz, sites)

    inter_hub = interaction_hub(frag, xyz, links)

    loop = loop_3d(frag, xyz)

    return None

if __name__ == "__main__":
    main() 
