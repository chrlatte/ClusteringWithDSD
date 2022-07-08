# """
# Author: Charlotte Versavel
# Date:   June 2022

#                             reworked_clusters_class.py

# Purpose: TODO -> WILL EVENTUALLY COPY IN FUNTIONS FROM CONNECTED_COMPONENTS_UTILS.PY -> FUNCTIONS CURRENTLY BEING DRAFTED IN CONNECTED_COMPONENTS_NOTEBOOK

# """

# from matrix_class import ProteinMatrix
# from matrix_class import SubMatrix
# from cluster_class import AllClusters
# from degreelist_class import DegreeList

# import numpy as np


# class ClusterInspection:

#     """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#     * * * * * * * * * * * * * MEMBER VARIABLES * * * * * * * * * * * * * *  
#     """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#     matrix: ProteinMatrix
#     clusters: AllClusters
#     degreelist: DegreeList

#     """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#     * * * * * * * * * * * * * * * FUNCTIONS * * * * * * * * * * * * * * *
#     """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

#     def __init__(self, protein_matrix: ProteinMatrix, all_clusters: AllClusters, degree_list: DegreeList) -> None:
#         """  
#         Parameters: TODO
#         Purpose:    TODO
#         Returns:    n/a
#         """
#         self.matrix = protein_matrix
#         self.clusters = all_clusters
#         self.degreelist = degree_list

    
#     def print_single_cluster_and_connected_proteins(self, x: int):
#         """
#         TODO
#         """
#         submatrix = SubMatrix(self.clusters.get_cluster_proteins(x), self.matrix)
#         num_components, labels = submatrix.get_num_components_and_labels()
        
#         print(f"Cluster {x} has {num_components} components: {[list(np.array(submatrix.get_list_of_proteins())[np.nonzero(labels == i)]) for i in range(num_components)]}.")

#         if num_components == 1: 
#             print(f"Cluster {x} is fully connected. will not search for attached proteins")
        
#         else:

#             component_dictionary = dict() # protein : component_num
#             j = 0
#             for array in [(np.array(submatrix.get_list_of_proteins())[np.nonzero(labels == i)]) for i in range(num_components)]:
#                 for protein in array:
#                     component_dictionary[protein] = j
#                 j += 1

#             list_of_proteins_connected_to_cluster = self.degreelist.create_list_of_proteins_connected_to_cluster(self.degreelist.get_list_of_proteins_sorted_by_degree(), self.clusters.get_cluster_proteins(x), min_num_connections=3)


#             for protein in list_of_proteins_connected_to_cluster:
#                 # print(f"{protein} is has 3+ connections to {matrix.get_list_of_proteins_connected_to_protein(protein)}")
#                 which_components = self.degreelist.which_components_of_a_cluster_would_a_protein_connect(protein, self.clusters.get_cluster_proteins(x), component_dictionary)

#                 if len(which_components) == num_components:
#                     print(f"protein {protein} has degree {self.matrix.find_degree(protein)} and will connect ALL {num_components} components!!!!!!!!!")
                
#                 elif len(which_components) > 1:
#                     print(f"protein {protein} has degree {self.matrix.find_degree(protein)} and will connect {len(which_components)} components: {which_components}")


#     def print_all_clusters_and_connected_proteins(self):
#         """
#         TODO
#         """
#         for x in range(self.clusters.get_num_clusters()):
#             print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
#             self.print_single_cluster_and_connected_proteins(x)


#     def find_clusters_that_match_criteria(self, ratio: float = 1/2, b: int = 0) -> list():
#         """
#         TODO
#         num_components must be less than ratio*num_proteins + b 
#         with a smaller ratio, or a lower b, the qualifying clusters will be more connected 


#         TODO TODO TODO -> copy in code from connected_components utils
#         """
#         pass
#         # clusters_that_match_criteria = list() # list of nums

#         # for x in range(self.clusters.get_num_clusters()):
#         #     submatrix = SubMatrix(self.clusters.get_cluster_proteins(x), self.matrix)

#         #     num_proteins = len(submatrix.get_list_of_proteins())
#         #     num_components, labels = submatrix.get_num_components_and_labels()
            
#         #     if (ratio * num_proteins + b) > num_components:
#         #         clusters_that_match_criteria.append(x)

        
#         # return clusters_that_match_criteria # TODO: return multiple things to pass into find_proteins function
        


#     def find_proteins_that_match_criteria(self, cluster_num: int, ratio: float = 1/2, b: int = 0, max_degree: int = 500) -> list():
#         """
#         TODO
#         a protein must connect more than ratio*num_components + b. if 

#         TODO : function could be improved by passing info from the find_cluster_that_match function, but for now, this is ok. 
#         """

#         submatrix = SubMatrix(self.clusters.get_cluster_proteins(cluster_num), self.matrix)
#         num_components, labels = submatrix.get_num_components_and_labels()


#         ### POPULATE COMPONENT DICTIONARY ###
#         component_dictionary = dict() # protein : component_num
#         j = 0
#         for array in [(np.array(submatrix.get_list_of_proteins())[np.nonzero(labels == i)]) for i in range(num_components)]:
#             for protein in array:
#                 component_dictionary[protein] = j
#             j += 1
#         # print(component_dictionary)

#         ## FIND CONNECTED PROTEINS AND DETERMINE IF THEY QUALIFY ###
#         qualifying_proteins = []
#         min_edges = 3

#         # print(f"self.degreelist.get_list_of_proteins_sorted_by_degree(): {self.degreelist.get_list_of_proteins_sorted_by_degree()}")

#         for protein in (self.degreelist.get_list_of_proteins_sorted_by_degree()):   
#             # print(f"in function! with protein {protein}")
#             print(f"protein {protein}. degree {self.matrix.find_degree(protein)}")

#             degree = self.matrix.find_degree(protein)
#             if (degree >= min_edges) and (degree <= max_degree):
#                 # print(f"after comparison!!. gonna determine num_edges: ")
                
#                 # print(f"calling fxn on line 146 with inputs protein: {protein}{type(protein)}, cluster_proteins: {list(self.clusters.get_cluster_proteins(cluster_num))} {type(list(self.clusters.get_cluster_proteins(cluster_num)))}, min_edges: {min_edges} {type(min_edges)}")
#                 num_edges, which_proteins = self.degreelist.determine_num_edges_to_cluster(protein, self.clusters.get_cluster_proteins(cluster_num), max_edges_until_return=min_edges, also_return_which_proteins=True)

#                 print(f"num edges : {num_edges}")

#                 if (num_edges >= min_edges): # if there are enough connections to a cluster, see how many different components a protein will connect
#                     print(f"PAST IF STATEMENT")
#                     which_components = self.degreelist.which_components_of_a_cluster_would_a_protein_connect(protein, self.clusters.get_cluster_proteins(cluster_num), component_dictionary, connected_proteins_within_cluster=which_proteins)

#                     if len(which_components) > ratio*num_components + b:
#                         print(f"SUCCESSSS: {protein} passes criteria for cluster {cluster_num}")
            


        


#         pass

#     def create_dict_of_proteins_to_add_to_clusters(self, x: int) -> dict():
#         """
#         returns a dict in a form  TODO
#         """
#         proteins_to_add_back = dict()

        

#         return proteins_to_add_back

