"""
Author: Charlotte Versavel
Date:   June 2022

                            reworked_clusters_class.py

Purpose: TODO

"""

from matrix_class import ProteinMatrix
from matrix_class import SubMatrix
from cluster_class import AllClusters
from degreelist_class import DegreeList

import numpy as np


class ClusterInspection:

    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    * * * * * * * * * * * * * MEMBER VARIABLES * * * * * * * * * * * * * *  
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    matrix: ProteinMatrix
    clusters: AllClusters
    degreelist: DegreeList

    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    * * * * * * * * * * * * * * * FUNCTIONS * * * * * * * * * * * * * * *
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    def __init__(self, protein_matrix: ProteinMatrix, all_clusters: AllClusters, degree_list: DegreeList) -> None:
        """  
        Parameters: TODO
        Purpose:    TODO
        Returns:    n/a
        """
        self.matrix = protein_matrix
        self.clusters = all_clusters
        self.degreelist = degree_list

    
    def print_single_cluster_and_connected_proteins(self, x: int):
        """
        TODO
        """
        submatrix = SubMatrix(self.clusters.get_cluster_proteins(x), self.matrix)
        num_components, labels = submatrix.get_num_components_and_labels()
        
        print(f"Cluster {x} has {num_components} components: {[list(np.array(submatrix.get_list_of_proteins())[np.nonzero(labels == i)]) for i in range(num_components)]}.")

        if num_components == 1: 
            print(f"Cluster {x} is fully connected. will not search for attached proteins")
        
        else:

            component_dictionary = dict() # protein : component_num
            j = 0
            for array in [(np.array(submatrix.get_list_of_proteins())[np.nonzero(labels == i)]) for i in range(num_components)]:
                for protein in array:
                    component_dictionary[protein] = j
                j += 1

            list_of_proteins_connected_to_cluster = self.degreelist.create_list_of_proteins_connected_to_cluster(self.degreelist.get_list_of_proteins_sorted_by_degree(), self.clusters.get_cluster_proteins(x), min_num_connections=3)


            for protein in list_of_proteins_connected_to_cluster:
                # print(f"{protein} is has 3+ connections to {matrix.get_list_of_proteins_connected_to_protein(protein)}")
                which_components = self.degreelist.which_components_of_a_cluster_would_a_protein_connect(protein, self.clusters.get_cluster_proteins(x), component_dictionary)

                if len(which_components) == num_components:
                    print(f"protein {protein} has degree {self.matrix.find_degree(protein)} and will connect ALL {num_components} components!!!!!!!!!")
                
                elif len(which_components) > 1:
                    print(f"protein {protein} has degree {self.matrix.find_degree(protein)} and will connect {len(which_components)} components: {which_components}")


    def print_all_clusters_and_connected_proteins(self):
        """
        TODO
        """
        for x in range(self.clusters.get_num_clusters()):
            print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
            self.print_single_cluster_and_connected_proteins(x)


    def find_clusters_that_match_criteria(self, ratio: float = 1/2, b: int = 0) -> list():
        """
        TODO
        num_components must be less than ratio*num_proteins + b 
        with a smaller ratio, or a lower b, the qualifying clusters will be more connected
        """
        clusters_that_match_criteria = list # list of nums

        for x in range(self.clusters.get_num_clusters()):

            submatrix = SubMatrix(self.clusters.get_cluster_proteins(x), self.matrix)
            num_components, labels = submatrix.get_num_components_and_labels()
            num_proteins = len(submatrix.get_list_of_proteins())

            

            if (ratio*num_proteins + b) > num_components:
                # print(f"cluster {x} meets criteria!! {num_components} components and {num_proteins} proteins")
                clusters_that_match_criteria.append(x)

            # else:
            #     print(f"cluster {x} has {num_components} components and {num_proteins} proteins")
        
        print(clusters_that_match_criteria)
        return clusters_that_match_criteria
        


    def find_proteins_that_match_criteria(self, x: int) -> list():
        """
        TODO
        """
        pass

    def create_dict_of_proteins_to_add_to_clusters(self, x: int) -> dict():
        """
        returns a dict in a form  TODO
        """
        proteins_to_add_back = dict()

        submatrix = SubMatrix(self.clusters.get_cluster_proteins(x), self.matrix)
        num_components, labels = submatrix.get_num_components_and_labels()
        
        print(f"Cluster {x} has {num_components} components: {[list(np.array(submatrix.get_list_of_proteins())[np.nonzero(labels == i)]) for i in range(num_components)]}.")

        if num_components == 1: 
            print(f"Cluster {x} is fully connected. will not search for attached proteins")
        
        else:

            component_dictionary = dict() # protein : component_num
            j = 0
            for array in [(np.array(submatrix.get_list_of_proteins())[np.nonzero(labels == i)]) for i in range(num_components)]:
                for protein in array:
                    component_dictionary[protein] = j
                j += 1

            list_of_proteins_connected_to_cluster = self.degreelist.create_list_of_proteins_connected_to_cluster(self.degreelist.get_list_of_proteins_sorted_by_degree(), self.clusters.get_cluster_proteins(x), min_num_connections=3)


            for protein in list_of_proteins_connected_to_cluster:
                # print(f"{protein} is has 3+ connections to {matrix.get_list_of_proteins_connected_to_protein(protein)}")
                which_components = self.degreelist.which_components_of_a_cluster_would_a_protein_connect(protein, self.clusters.get_cluster_proteins(x), component_dictionary)

                if len(which_components) == num_components:
                    print(f"protein {protein} has degree {self.matrix.find_degree(protein)} and will connect ALL {num_components} components!!!!!!!!!")
                
                elif len(which_components) > 1:
                    print(f"protein {protein} has degree {self.matrix.find_degree(protein)} and will connect {len(which_components)} components: {which_components}")



        return proteins_to_add_back

