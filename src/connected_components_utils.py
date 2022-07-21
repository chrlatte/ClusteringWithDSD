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





def find_clusters_and_proteins_together():
    """
    function is faster
    """