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
    # print(f"Matrix:\n{matrix}")


    clusters = ProteinClusters("../data/testing_data/fake_cluster.txt")
    # print(f"Clusters:\n{clusters}")
    # clusters.print_all()


    degreelist = DegreeList(matrix)
    # print(f"Degree list:\n{degreelist}")

    # print(f"cluster 2:\n{clusters.get_cluster(2)}")

    # matrix.has_edge("PRKCA", "KRAS")

    # degreelist.determine_num_edges_to_cluster("PRKCA", ["KRAS", "CDKN1A"])

    #print(f"list of proteins: {matrix.get_list_of_proteins}")
    print(matrix.get_list_of_proteins)
    #degreelist.create_list_of_proteins_connected_to_cluster(matrix.get_list_of_proteins, clusters.get_proteins_from_cluster(0))
    

    

     
    


if __name__ == "__main__":
    main()
