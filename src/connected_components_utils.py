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
        return .95
    elif num_clusters > 50:
        return .995
    else: 
        return 1
    