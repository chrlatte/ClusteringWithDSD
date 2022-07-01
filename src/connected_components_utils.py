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
def print_single_cluster_and_connected_proteins(x: int, matrix: ProteinMatrix, clusters: AllClusters, degreelist: DegreeList):
    """
    TODO
    """
    submatrix = SubMatrix(clusters.get_cluster_proteins(x), matrix)
    num_components, labels = submatrix.get_num_components_and_labels()
    
    print(f"Cluster {x} has {num_components} components: {[list(np.array(submatrix.get_list_of_proteins())[np.nonzero(labels == i)]) for i in range(num_components)]}.")

    if num_components == 1: 
        print(f"Cluster {x} is fully connected. will not search for attached proteins")
    
    else:

        component_dictionary = dict() # protein : component_num
        j = 0
        for array in [(np.array(submatrix.get_list_of_proteins())[np.nonzero(labels == i)]) for i in range(num_components)]:
            for protein in array:
                component_dictionary[protein] = j
            j += 1

        list_of_proteins_connected_to_cluster = degreelist.create_list_of_proteins_connected_to_cluster(degreelist.get_list_of_proteins_sorted_by_degree(), clusters.get_cluster_proteins(x), min_num_connections=3)


        for protein in list_of_proteins_connected_to_cluster:
            # print(f"{protein} is has 3+ connections to {matrix.get_list_of_proteins_connected_to_protein(protein)}")
            which_components = degreelist.which_components_of_a_cluster_would_a_protein_connect(protein, clusters.get_cluster_proteins(x), component_dictionary)

            if len(which_components) == num_components:
                print(f"protein {protein} has degree {matrix.find_degree(protein)} and will connect ALL {num_components} components!!!!!!!!!")
            
            elif len(which_components) > 1:
                print(f"protein {protein} has degree {matrix.find_degree(protein)} and will connect {len(which_components)} components: {which_components}")


def print_all_clusters_and_connected_proteins(matrix: ProteinMatrix, clusters: AllClusters, degreelist: DegreeList):
    """
    TODO
    """
    for x in range(clusters.get_num_clusters()):
        print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
        print_single_cluster_and_connected_proteins(x, matrix, clusters, degreelist)


def find_clusters_that_match_criteria(matrix: ProteinMatrix, clusters: AllClusters, degreelist: DegreeList) -> list():
    """
    TODO
    """
    pass

def find_proteins_that_match_criteria(x: int, matrix: ProteinMatrix, clusters: AllClusters, degreelist: DegreeList) -> list():
    """
    TODO
    """
    pass

def create_dict_of_proteins_to_add_to_clusters(x: int, matrix: ProteinMatrix, clusters: AllClusters, degreelist: DegreeList) -> dict():
    """
    returns a dict in a form  TODO
    """
    proteins_to_add_back = dict()

    submatrix = SubMatrix(clusters.get_cluster_proteins(x), matrix)
    num_components, labels = submatrix.get_num_components_and_labels()
    
    print(f"Cluster {x} has {num_components} components: {[list(np.array(submatrix.get_list_of_proteins())[np.nonzero(labels == i)]) for i in range(num_components)]}.")

    if num_components == 1: 
        print(f"Cluster {x} is fully connected. will not search for attached proteins")
    
    else:

        component_dictionary = dict() # protein : component_num
        j = 0
        for array in [(np.array(submatrix.get_list_of_proteins())[np.nonzero(labels == i)]) for i in range(num_components)]:
            for protein in array:
                component_dictionary[protein] = j
            j += 1

        list_of_proteins_connected_to_cluster = degreelist.create_list_of_proteins_connected_to_cluster(degreelist.get_list_of_proteins_sorted_by_degree(), clusters.get_cluster_proteins(x), min_num_connections=3)


        for protein in list_of_proteins_connected_to_cluster:
            # print(f"{protein} is has 3+ connections to {matrix.get_list_of_proteins_connected_to_protein(protein)}")
            which_components = degreelist.which_components_of_a_cluster_would_a_protein_connect(protein, clusters.get_cluster_proteins(x), component_dictionary)

            if len(which_components) == num_components:
                print(f"protein {protein} has degree {matrix.find_degree(protein)} and will connect ALL {num_components} components!!!!!!!!!")
            
            elif len(which_components) > 1:
                print(f"protein {protein} has degree {matrix.find_degree(protein)} and will connect {len(which_components)} components: {which_components}")



    return proteins_to_add_back

