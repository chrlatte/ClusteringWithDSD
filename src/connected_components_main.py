"""
Author: Charlotte Versavel
Date:   June 2022

                            connected_components.py

Purpose: a main to use the scipy library to determine connected components in a cluster

"""

from matrix_class import *
from cluster_class import *
from degreelist_class import *
import sys

# from scipy.sparse import coo_matrix
# from scipy.sparse import csr_matrix
# from scipy.sparse.csgraph import connected_components


def main():

    smaller_testing_matrix_file = "../data/testing_data/tiny_dream3.txt"
    testing_matrix_file = "../data/testing_data/small_dream3.txt"

    matrix_file_for_cluster = "../data/testing_data/fake_cluster_dream.txt"
    testing_cluster_file = "../data/testing_data/fake_cluster.txt"

    actual_matrix_file = "../data/networks/DREAM_files/dream_2.txt"
    smaller_cluster_file = "../data/testing_data/moderately_connected_clusters.txt"
    actual_cluster_file = "../data/clusters/3344522.7320912.1_ppi_anonym_v2.txt"


    matrix = ProteinMatrix(matrix_file_for_cluster)
    # print(f"Matrix:\n{matrix}")

    clusters = AllClusters(testing_cluster_file)
    # print(f"Clusters:\n{clusters}")
    clusters.print_all()

    degreelist = DegreeList(matrix)
    # print(f"Degree list:\n{degreelist}")

 
    

    for i in range(clusters.get_num_clusters()):

        print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")

        submatrix = SubMatrix(clusters.get_cluster_proteins(i), matrix)
        # print(f"SUBMATRIX FROM CLUSTER {i}: \n{submatrix.get_matrix()}")
        n, labels = submatrix.get_num_components_and_labels()
        print(f"Cluster {i} has {n} components: {[list(np.array(submatrix.get_list_of_proteins())[np.nonzero(labels == i)]) for i in range(n)]}.")
        
        list_of_proteins_connected_to_cluster = degreelist.create_list_of_proteins_connected_to_cluster(degreelist.get_list_of_proteins_sorted_by_degree(), clusters.get_cluster_proteins(i), min_num_connections=3)

        print(f"proteins connected 3+ times to cluster {i}: {list_of_proteins_connected_to_cluster}")

        
        for protein in list_of_proteins_connected_to_cluster:

            will_connect = degreelist.determine_if_a_protein_will_connect_a_cluster(protein, clusters.get_cluster_proteins(i), labels)
            
            if will_connect:
                print(f"protein {protein} will connect components {degreelist.which_components_of_a_cluster_would_a_protein_connect(protein, clusters.get_cluster_proteins(i), labels)}")
            else:
                print(f"{protein} will not connect cluster.")



if __name__ == "__main__":
    main()
