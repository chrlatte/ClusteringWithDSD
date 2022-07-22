"""
Author: Charlotte Versavel
Date:   July 2022

                            analysis_main.py

Purpose: a main TODO

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


    # Dream3 Data:
    clusters = AllClusters(protein_to_cluster_dict=dream3_clusters_dict)
    matrix = ProteinMatrix(dream3_matrix_file)
    degreelist = DegreeList(matrix)

    # print(f"Matrix:\n{matrix}")
    # print(f"Clusters:\n{clusters}")
    # clusters.print_all()
    # print(f"Degree list:\n{degreelist}")


    qualifying_clusters, qualifying_proteins = find_clusters_and_proteins_together(matrix, clusters, degreelist, cluster_ratio=pick_ratio(clusters.get_num_clusters()), protein_ratio=.05, protein_constant=2)

    print(f"{qualifying_proteins}")
    
    

if __name__ == "__main__":
    main()
