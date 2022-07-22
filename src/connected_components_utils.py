"""
Author: Charlotte Versavel
Date:   June 2022

                            connected_components_utils.py

Purpose: TODO

"""

from matrix_class import ProteinMatrix
from matrix_class import SubMatrix
from cluster_class import AllClusters
from degreelist_class import DegreeList

import numpy as np

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
 * * * * * * * * * * * * * * * FUNCTIONS * * * * * * * * * * * * * * *
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
def find_clusters_that_match_criteria(matrix: ProteinMatrix, clusters: AllClusters, degreelist: DegreeList, ratio: float = .5, constant: int = 0) -> list():
    """
    Parameters: 
        matrix - a ProteinMatrix of all protein interactions
        clusters - an AllClusters containing proteins grouped into clusters
        ratio and constant - used together to determine which clusters qualify, with the output of the function being ratio * input + constant
    Purpose:    function goes through all clusters and finds those that are 
                'connected enough'. to qualify, the number of components in the 
                cluster must be less the value of the qualifying_threshhold 
                function when the number of proteins in the cluster is passed in
    Returns:    a list containing the numbers of the clusters that qualify
    """
    
    cluster_nums_that_qualify = list()

    for key in clusters.get_all_clusters():
        # create a submatrix out of the proteins in the cluster
        submatrix = SubMatrix(clusters.get_cluster_proteins(key), matrix)
        num_components, labels = submatrix.get_num_components_and_labels()
        # print(f"num components is {num_components}. num proteins is {len(submatrix.get_list_of_proteins())}")
        if num_components < ratio * len(submatrix.get_list_of_proteins()) + constant:
            # print('success')
            cluster_nums_that_qualify.append(key)

    return cluster_nums_that_qualify


def pick_ratio(num_clusters: int):
    """
    will determine an approximate ratio to start with based on the total number of clusters
    """
    if (num_clusters > 1000):
        return .5
    elif num_clusters > 500:
        return .7
    elif num_clusters > 200:
        return .9
    elif num_clusters > 100:
        return .925
    elif num_clusters > 50:
        return .995
    else: 
        return 1


def find_proteins_that_match_criteria(cluster_num: int, matrix: ProteinMatrix, clusters: AllClusters, degreelist: DegreeList, ratio: float = .5, constant: int = 0, min_components_that_protein_connects: int = -1, max_degree: int = 500) -> list():
    """
    can choose to a ratio and a constant, or a min num components
    a protein must connect more than ratio*num_components + constant
    TODO : function could be improved by passing (submatrix) info from the find_cluster_that_match function, but for now, this is ok. (so that it doesn't need to be reconstructed)
    """
    if (min_components_that_protein_connects == -1):
            min_components_that_protein_connects = constant + ratio * len(clusters.get_cluster_proteins(cluster_num))
        
    submatrix = SubMatrix(clusters.get_cluster_proteins(cluster_num), matrix)
    num_components, labels = submatrix.get_num_components_and_labels()


    ### POPULATE COMPONENT DICTIONARY ###
    component_dictionary = dict() # protein : component_num
    j = 0
    for array in [(np.array(submatrix.get_list_of_proteins())[np.nonzero(labels == i)]) for i in range(num_components)]:
        for protein in array:
            component_dictionary[protein] = j
        j += 1
    
    ## FIND CONNECTED PROTEINS AND DETERMINE IF THEY QUALIFY ###
    qualifying_proteins = list()

    for protein in (degreelist.get_list_of_proteins_sorted_by_degree()):   
        degree = matrix.find_degree(protein)

        if (degree >= min_components_that_protein_connects) and (degree <= max_degree):
            num_edges, which_proteins = degreelist.determine_num_edges_to_cluster(protein, clusters.get_cluster_proteins(cluster_num), also_return_which_proteins=True)
                
            if (num_edges >= min_components_that_protein_connects):
                set_of_components_that_protein_connects = degreelist.which_components_of_a_cluster_would_a_protein_connect(protein, clusters.get_cluster_proteins(cluster_num), component_dictionary, connected_proteins_within_cluster=which_proteins)

                if len(set_of_components_that_protein_connects) >= min_components_that_protein_connects:
                    qualifying_proteins.append(protein)

    return qualifying_proteins


def qualifying_proteins_using_submatrix(cluster_num: int, submatrix: SubMatrix, clusters: AllClusters, degreelist: DegreeList, ratio: float = .5, constant: int = 0, min_components_that_protein_connects: int = -1, max_degree: int = 500) -> list():
    """
    TODO : a revised version of the find_proteins_that_match_criteria function that takes in a submatrix as a parameter, and therefore doesn't need to construct one. 
    """
    if (min_components_that_protein_connects == -1):
            min_components_that_protein_connects = constant + ratio * len(clusters.get_cluster_proteins(cluster_num))
        
    num_components, labels = submatrix.get_num_components_and_labels()

    ### POPULATE COMPONENT DICTIONARY ###
    component_dictionary = dict() # protein : component_num
    j = 0
    for array in [(np.array(submatrix.get_list_of_proteins())[np.nonzero(labels == i)]) for i in range(num_components)]:
        for protein in array:
            component_dictionary[protein] = j
        j += 1
    
    ## FIND CONNECTED PROTEINS AND DETERMINE IF THEY QUALIFY 
    qualifying_proteins = list()

    for protein in (degreelist.get_list_of_proteins_sorted_by_degree()):   
        num_edges, which_proteins = degreelist.determine_num_edges_to_cluster(protein, clusters.get_cluster_proteins(cluster_num), also_return_which_proteins=True)
                
        if (num_edges >= min_components_that_protein_connects):
            set_of_components_that_protein_connects = degreelist.which_components_of_a_cluster_would_a_protein_connect(protein, clusters.get_cluster_proteins(cluster_num), component_dictionary, connected_proteins_within_cluster=which_proteins)

            if len(set_of_components_that_protein_connects) >= min_components_that_protein_connects:
                qualifying_proteins.append(protein)

    return qualifying_proteins



def find_clusters_and_proteins_together(matrix: ProteinMatrix, clusters: AllClusters, degreelist: DegreeList, cluster_ratio: float = .5, cluster_constant: int = 0, protein_ratio: float = .5, protein_constant: int = 0, min_components_that_protein_connects: int = -1, max_degree: int = 500) -> list() and dict():
    """
    function is a version of find_clusters_that_match_criteria, that, once it finds the cluster, finds corresponding proteins at the same time so that the submatrix doesn't need to be reconstructed

    Parameters: 
        matrix - a ProteinMatrix of all protein interactions
        clusters - an AllClusters containing proteins grouped into clusters
        cluster_ratio and cluster_constant - used together to determine which clusters qualify, with the output of the function being cluster_ratio * input + cluster_constant
        TODO: remaining parameters
    Purpose:    determines clusters that are mostly highly connected, then 
                determines which proteins that, when added to the cluster, will 
                increase it's connectedness
    Returns:    a list containing the numbers of the clusters that qualify, and 
                a dictionary linking each cluster, to a list of the qualifying 
                proteins
    """
    
    cluster_nums_that_qualify = list()
    qualifying_proteins_dict = dict()

    for cluster_num in clusters.get_all_clusters():
        # create a submatrix out of the proteins in the cluster
        submatrix = SubMatrix(clusters.get_cluster_proteins(cluster_num), matrix)
        num_components, labels = submatrix.get_num_components_and_labels()
        # print(f"num components is {num_components}. num proteins is {len(submatrix.get_list_of_proteins())}")
        if num_components < cluster_ratio * len(submatrix.get_list_of_proteins()) + cluster_constant:

            # add cluster to list showing that it qualifies, 
            cluster_nums_that_qualify.append(cluster_num)
            # then do analysis on the cluster
            qualifying_proteins_dict[cluster_num] = qualifying_proteins_using_submatrix(cluster_num, submatrix, clusters, degreelist, ratio=protein_ratio, constant=protein_constant, min_components_that_protein_connects=min_components_that_protein_connects, max_degree=max_degree)


    return cluster_nums_that_qualify, qualifying_proteins_dict


def create_new_clusters(clusters_to_qualifying_proteins: dict(), csv_filename: str = "", protein_to_cluster_dict: dict() = {}, original_clusters: AllClusters = AllClusters()) -> AllClusters:
    """
    csv_filename: str = "", protein_to_cluster_dict: dict() ={}, original_clusters: AllClusters = AllClusters(), are all different ways to pass in info to make new clusters (and you should choose one of them). Please note, that if you use Original Clusters, the original clusters will be modified to include the new qualifying proteins
    """
    modified_clusters = AllClusters()
    
    if csv_filename != "":
        modified_clusters = AllClusters(csv_filename=csv_filename)
    elif protein_to_cluster_dict: # dictionary not empty
        modified_clusters = AllClusters(protein_to_cluster_dict=protein_to_cluster_dict)
    else:
        modified_clusters = original_clusters
    
    for key in clusters_to_qualifying_proteins:
        if clusters_to_qualifying_proteins[key]: # will return true if this cluster has qualifying proteins (list not empty)
            for protein in clusters_to_qualifying_proteins[key]:
                modified_clusters.add_protein_to_cluster(protein, key)
    return modified_clusters