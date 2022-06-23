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
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components


def main():

    testing_matrix_file = "../data/testing_data/small_dream3.txt"
    testing_cluster_file = "../data/testing_data/fake_cluster.txt"

    smaller_testing_matrix_file = "../data/testing_data/tiny_dream3.txt"


    actual_matrix_file = "../data/networks/DREAM_files/dream_3.txt"
    actual_cluster_file = "../data/clusters/3344522.7320912.1_ppi_anonym_v2.txt"

    matrix = ProteinMatrix(smaller_testing_matrix_file)
    # print(f"Matrix:\n{matrix}")

    clusters = ProteinClusters(testing_cluster_file)
    # print(f"Clusters:\n{clusters}")
    # clusters.print_all()

    degreelist = DegreeList(matrix)
    # print(f"Degree list:\n{degreelist}")

    
    # graph = matrix.get_matrix()
    # print(graph)

    # want to take a cluster and then make a submatrix. maybe the submatrix can 
    # have the csr matrix
    clusters.get_cluster(0)



    # from pandas.api.types import CategoricalDtype
    # from scipy import sparse

    # sparse_df = matrix.get_matrix().astype(pd.SparseDtype(object, fill_value=0))

    # proteins_cat = CategoricalDtype(categories=sorted(matrix.get_list_of_proteins()), ordered=True)

    # movies = df["movie_id"].unique()
    # shape = (len(users), len(movies))

    # # Create indices for users and movies
    # user_cat = CategoricalDtype(categories=sorted(users), ordered=True)
    # movie_cat = CategoricalDtype(categories=sorted(movies), ordered=True)
    # user_index = df["user_id"].astype(user_cat).cat.codes
    # movie_index = df["movie_id"].astype(movie_cat).cat.codes

    # print(matrix.get_matrix())
    # matrix2 = coo_matrix(matrix.get_matrix())
    # print(matrix2)
    # print(type(matrix2))
    # csr = matrix2.tocsr()
    # print(csr)


    # # Conversion via COO matrix
    # coo = sparse.coo_matrix((df["rating"], (user_index, movie_index)), shape=shape)
    # csr = coo.tocsr()






    

    matrix_csr = csr_matrix(matrix.get_matrix().astype(pd.SparseDtype(dtype=float, fill_value=0)))
    print(matrix_csr)

    n_components, labels = connected_components(matrix_csr, directed=False, return_labels=True)
    print(f"n_components: {n_components}; labels: {labels}")

    # # print(df_sparsed)
    # csr_matrix(df_sparsed)
    # # csr_matrix(df_sparsed) # Note you need .sparse accessor to access .to_coo()

    # graph = csr_matrix(graph)
    # print(graph)
    #n_components, labels = connected_components(csgraph=graph, directed=False, return_labels=True)
    
    
    # print(degreelist.get_list_of_proteins_sorted_by_degree())
    #########
    
    # individual_cluster = clusters.get_cluster(0)
    # clusters_matrix = SubMatrix(individual_cluster.get_protein_list(), matrix)
    # print()


    
    # for n in range(clusters.get_num_clusters()):
    #     result = degreelist.create_list_of_proteins_connected_to_cluster(degreelist.get_list_of_proteins_sorted_by_degree(), clusters.get_proteins_from_cluster(n))
    #     print(f"THE PROTEINS WITH 3+ CONNECTIONS TO CLUSTER {n} ARE {result}")



    # print(f"using cluster: {individual_cluster}")
    # print(f"using list of proteins: {matrix.get_list_of_proteins()}")
    # print(f"visual representation of cluster:\n{clusters_matrix.get_matrix()}")

    # print(f"THE PROTEINS WITH >3 CONNECTIONS TO CLUSTER: {individual_cluster} ARE {result}")
    

    

     
    


if __name__ == "__main__":
    main()
