"""
Author: Charlotte Versavel
Date:   July 2022

                            reworked_clusters_main.py

Purpose: TODO

"""

import numpy as np
import json


from matrix_class import ProteinMatrix
from matrix_class import SubMatrix
from cluster_class import AllClusters
from degreelist_class import DegreeList

# from reworked_clusters_class import ClusterInspection

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
    # matrix = ProteinMatrix(testing_matrix_file)
    # print(f"Matrix:\n{matrix}")


    
    clusters = AllClusters(protein_to_cluster_dict=dict_of_clusters)
    # clusters = AllClusters(testing_cluster_file)


    # # print(f"Clusters:\n{clusters}")
    # clusters.print_all()

    degreelist = DegreeList(matrix)
    # print(f"Degree list:\n{degreelist}")

    clusters_that_are_somewhat_connected = find_clusters_that_match_criteria(matrix, clusters, degreelist, ratio=1)
    print(f"clusters_that_are_somewhat_connected: {clusters_that_are_somewhat_connected}")


    
if __name__ == "__main__":
    main()
