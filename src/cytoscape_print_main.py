"""
Author: Charlotte Versavel
Date:   June 2022

                          cytoscape_print_main.py

Purpose: a main to print to a .sif file that can be uploaded to cytoscape

"""

from cluster_class import *
from degreelist_class import *
from cytoscape_print_class import *


def main():

    testing_matrix_file = "../data/testing_data/fake_cluster_dream.txt"
    testing_cluster_file = "../data/testing_data/fake_cluster.txt"

    dream3_matrix_file = "../data/networks/DREAM_files/dream_3.txt"
    dream3_old_cluster_file = "../data/clusters/3344522.7320912.1_ppi_anonym_v2.txt"


    
    matrix = ProteinMatrix(dream3_matrix_file)
    clusters = AllClusters(dream3_old_cluster_file)

    foo = PrintToFile()

    
    foo.print_all_interactions(matrix)
    foo.assign_colors_to_clusters(clusters)





    print(f"done")

if __name__ == "__main__":
    main()
