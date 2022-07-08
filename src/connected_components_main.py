"""
Author: Charlotte Versavel
Date:   June 2022

                            connected_components.py

Purpose: a main to use the scipy library to determine connected components in a cluster

"""

import numpy as np
import json


from matrix_class import ProteinMatrix
from matrix_class import SubMatrix
from cluster_class import AllClusters
from degreelist_class import DegreeList


from connected_components_utils import *


def main():

    testing_matrix_file = "../data/testing_data/fake_cluster_dream.txt"
    testing_cluster_file = "../data/testing_data/fake_cluster.txt"

    dream3_matrix_file = "../data/networks/DREAM_files/dream_3.txt"
    
    dream3_cluster_dict = "../data/results/DREAM-3-cc/d3_5_100.json-cluster.json"
    dict_of_clusters = {}
    # # convert actual cluster file to a dictionary!!
    with open(dream3_cluster_dict,"r") as cluster_dict_file:
        dict_of_clusters = json.load(cluster_dict_file)
    

    dream2_matrix_file = "../data/networks/DREAM_files/dream_2.txt"
    dream3_old_cluster_file = "../data/clusters/3344522.7320912.1_ppi_anonym_v2.txt"


    matrix = ProteinMatrix(dream3_matrix_file)
    # print(f"Matrix:\n{matrix}")


    
    clusters = AllClusters(protein_to_cluster_dict=dict_of_clusters)
    # # print(f"Clusters:\n{clusters}")
    # clusters.print_all()

    degreelist = DegreeList(matrix)
    # print(f"Degree list:\n{degreelist}")

 
    

    # print_cluster_components_and_proteins_that_are_connected(matrix, clusters, degreelist)
    # clusters.get_cluster_proteins(0) # TODO is empty

    
    # print_all_clusters_and_connected_proteins(matrix, clusters, degreelist)
    # create_dict_of_proteins_to_add_to_clusters()
    clusters_that_are_somewhat_connected = find_clusters_that_match_criteria(matrix, clusters, degreelist, ratio=1)
    print(f"clusters_that_are_somewhat_connected: {clusters_that_are_somewhat_connected}")

if __name__ == "__main__":
    main()
