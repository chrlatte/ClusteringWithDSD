"""
Author: Charlotte Versavel
Date:   June 2022

                            connected_components.py

Purpose: a main to use the scipy library to determine connected components in a cluster

"""

from matrix_class import *
from cluster_class import *
from degreelist_class import *

# from scipy.sparse import coo_matrix
# from scipy.sparse import csr_matrix
# from scipy.sparse.csgraph import connected_components


def main():

    smaller_testing_matrix_file = "../data/testing_data/tiny_dream3.txt"
    testing_matrix_file = "../data/testing_data/small_dream3.txt"
    testing_cluster_file = "../data/testing_data/fake_cluster.txt"

    actual_matrix_file = "../data/networks/DREAM_files/dream_2.txt"
    smaller_cluster_file = "../data/testing_data/moderately_connected_clusters.txt"
    actual_cluster_file = "../data/clusters/3344522.7320912.1_ppi_anonym_v2.txt"


    matrix = ProteinMatrix(actual_matrix_file)
    # print(f"Matrix:\n{matrix}")

    clusters = ProteinClusters(actual_cluster_file)
    # print(f"Clusters:\n{clusters}")
    # clusters.print_all()

    degreelist = DegreeList(matrix)
    # print(f"Degree list:\n{degreelist}")


    # want to take a cluster and then make a submatrix. 
    # the submatrix houses the CSR matrix
    for i in range(clusters.get_num_clusters()): # clusters.get_num_clusters()

        submatrix = SubMatrix(clusters.get_cluster_proteins(i), matrix)
        # print(submatrix.get_matrix())
        n, labels = submatrix.get_num_components_and_labels()
        print(f"Cluster {i} has {n} components: {[list((submatrix.get_list_of_proteins())[np.nonzero(labels == i)]) for i in range(n)]}\n")



if __name__ == "__main__":
    main()
