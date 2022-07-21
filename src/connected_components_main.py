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
    dream3_cluster_file = "../data/results/DREAM-3-cc/d3_5_100.json-cluster.json" 

    dream3_clusters_dict = {}
    # convert actual cluster file to a dictionary!!
    with open(dream3_cluster_file,"r") as cluster_dict_file:
        dream3_clusters_dict = json.load(cluster_dict_file)
    

    # # Testing Data:
    # matrix = ProteinMatrix(testing_matrix_file)
    # clusters = AllClusters(testing_cluster_file)
    # degreelist = DegreeList(matrix)

    # Dream3 Data:
    clusters = AllClusters(protein_to_cluster_dict=dream3_clusters_dict)
    matrix = ProteinMatrix(dream3_matrix_file)
    degreelist = DegreeList(matrix)

    # print(f"Matrix:\n{matrix}")
    # print(f"Clusters:\n{clusters}")
    # clusters.print_all()
    # print(f"Degree list:\n{degreelist}")


    # first, going to find the clusters that are moderatly connected (by dream3 standards)
    clusters_that_are_somewhat_connected = find_clusters_that_match_criteria(matrix, clusters, degreelist, ratio=pick_ratio(clusters.get_num_clusters()))

    print(f"clusters_that_are_somewhat_connected: {clusters_that_are_somewhat_connected}")

    

if __name__ == "__main__":
    main()
