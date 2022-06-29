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
    # dream3_matrix_file = "../data/networks/DREAM_files/dream_3.txt"

    matrix = ProteinMatrix(testing_matrix_file)
    

    sif = PrintToFile()

    sif.clear_file()
    sif.print_matrix(matrix)

    sif.swap_file("cytoscape_colors.txt")
    sif.clear_file()

    sif.assign_colors_to_clusters(["C0P1", "C1P1", "C2P1"])
    sif.assign_colors_to_clusters(["C0P2", "C1P2", "C2P2"])




    print(f"done")

if __name__ == "__main__":
    main()
