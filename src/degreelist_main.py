"""
Author: Charlotte Versavel
Date:   June 2022

                             degreelist_main.py

Purpose: a main to test the matrix and cluster classes together, focusing mainly on the degreelist generation and implementation

"""

from matrix_class import *
from cluster_class import *
from degreelist_class import *

def main():

    matrix = ProteinMatrix("../data/testing_data/tiny_dream3.txt")
    print(f"Matrix:\n{matrix}")


    clusters = ProteinClusters("../data/testing_data/fake_cluster.txt")
    print(f"Clusters:\n{clusters}")
    clusters.print_all()


    degreelist = DegreeList(matrix)
    print(f"Degree list:\n{degreelist}")
    

    

     
    


if __name__ == "__main__":
    main()
