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

    matrix = ProteinMatrix("../data/networks/DREAM_files/dream_3.txt")
    # print(f"Matrix:\n{matrix}")

    clusters = ProteinClusters("../data/clusters/3344522.7320912.1_ppi_anonym_v2.txt")
    # print(f"Clusters:\n{clusters}")
    # clusters.print_all()


    degreelist = DegreeList(matrix)
    # print(f"Degree list:\n{degreelist}")


    #########
    
    # individual_cluster = clusters.get_cluster(0)
    # clusters_matrix = SubMatrix(individual_cluster.get_protein_list(), matrix)
    # print()


    
    for n in range(clusters.get_num_clusters()):
        result = degreelist.create_list_of_proteins_connected_to_cluster(matrix.get_list_of_proteins(), clusters.get_proteins_from_cluster(n))
        print(f"THE PROTEINS WITH >3 CONNECTIONS TO CLUSTER: {n} ARE {result}")



    # print(f"using cluster: {individual_cluster}")
    # print(f"using list of proteins: {matrix.get_list_of_proteins()}")
    # print(f"visual representation of cluster:\n{clusters_matrix.get_matrix()}")

    # print(f"THE PROTEINS WITH >3 CONNECTIONS TO CLUSTER: {individual_cluster} ARE {result}")
    

    

     
    


if __name__ == "__main__":
    main()
