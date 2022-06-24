"""
Author: Charlotte Versavel
Date:   June 2022

                             degreelist_class.py

Purpose: a class TODO

"""

import pandas as pd 
import numpy as np

from matrix_class import *



class DegreeList:
    
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    * * * * * * * * * * * * * MEMBER VARIABLES * * * * * * * * * * * * * *  
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    protein_matrix : ProteinMatrix = ProteinMatrix

    sorted_protein_degree_dict = dict()


    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    * * * * * * * * * * * * * * INITIALIZERS * * * * * * * * * * * * * * *  
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    def __init__(self, matrix : ProteinMatrix) -> None:
        """            
        Parameters: matrix is populated with proteins and their interaction 
                    weights
        Purpose:    to take in a proteinMatrix (or submatrix) and create a 
                    sorted dictionary of protein:degree for all proteins in the 
                    matrix.
        Returns:    n/a
        """
        self.protein_matrix = matrix

        protein_degree_dict = {name:matrix.find_degree(name) for name in matrix.get_list_of_proteins()}

        self.sorted_protein_degree_dict = sorted(protein_degree_dict.items(), key=lambda x: x[1], reverse=True)


    def __repr__(self): 
        """             
        Purpose:    to override the print function for this class to print the 
                    sorted dictionary when called
        Returns:    a string of the dictionary
        """
        
        return str(self.sorted_protein_degree_dict)

    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    * * * * * * * * * * * * * * * GETTERS * * * * * * * * * * * * * * * * *  
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    def get_degree_list(self) -> list():
        """             
        Purpose:    to allow access to the sorted degree list
        Returns:    a list of tuples of (protein, degree)
        """
        return self.sorted_protein_degree_dict
    
    def get_list_of_proteins_sorted_by_degree(self) -> list():
         
        reverse_list_of_proteins = []

        for protein_degree_pair in reversed(self.sorted_protein_degree_dict):
            reverse_list_of_proteins.append(protein_degree_pair[0])

        return reverse_list_of_proteins

    def get_protein_at_index(self, index : int, degree = False) -> str or tuple:
        """             
        Parameters: index is the index of the protein in the sorted list
                    degree is a boolean that determines if the degree is returned as well
        Purpose:    to return the protein at the specified index
        Returns:    the protein at the specified index, or if degree is True, a 
                    tuple of (protein, degree)
        """
        if not degree:
            return self.sorted_protein_degree_dict[index][0]
        else:
            return self.sorted_protein_degree_dict[index]
    

    
    
    def determine_num_edges_to_cluster(self, protein : str, cluster_list : np.array, max_edges_until_return : int = -1, also_return_which_proteins : bool = False) -> int:
        """             
        Parameters: protein is a single protein in the matrix
                    cluster_list is a list of proteins in a cluster
                    max_edges_until_return allows the function to stop counting edges once a certain target is reached
                    also_return_which_proteins if set to true, will return which proteins have have edges to the given protein. should not be set to true when max_edges is set to a specific value
        Purpose:    to determine the number of edges between the protein and the proteins in the cluster
        Returns:    the number of edges
        """
        for already_in_cluster in cluster_list:
            if protein == already_in_cluster:
                return 0
        
        num_edges = 0
        which_proteins = list() 
        # print(cluster_list)
        if max_edges_until_return == -1: # max_edges_until_return has been left unspecified
            for cluster_protein in cluster_list:

                if (self.protein_matrix).has_edge(protein, cluster_protein):
                    num_edges += 1
                    which_proteins.append(cluster_protein)
        else: # max_edges_until_return has been specified
            for cluster_protein in cluster_list:
                # print(f"about to call has edge 1")
                if (self.protein_matrix).has_edge(protein, cluster_protein):
                    num_edges += 1
                if num_edges >= max_edges_until_return:
                    return num_edges
        
        
        if (also_return_which_proteins):
            return num_edges, which_proteins
        return num_edges
        


    def create_list_of_proteins_connected_to_cluster(self, list_of_proteins: np.array, cluster_list : np.array, max_list_length : int = -1, min_num_connections : int = 3) -> list:
        """             
        Parameters: cluster_list is a list of proteins in a cluster
                    max_list_length is an upper bound for the number of proteins to return in a list. If None, all proteins with at least min_num_connections connections are added to the list
                    min_num_connections is the minimum number of connections a protein must have to be added to the list and considered 'connected' to the cluster
        Purpose:    to create a list of proteins that are connected to the cluster
        Returns:    a list of proteins that are connected to the cluster
        """
        
        qualifying_proteins = []
        if max_list_length == -1: # max_list_length left unspecified
            for protein in list_of_proteins:
                num_edges = self.determine_num_edges_to_cluster(protein, cluster_list, max_edges_until_return=min_num_connections)

                if (num_edges >= min_num_connections):
                        qualifying_proteins.append(protein)
        else: # max_list_length has been specified    
            for protein in list_of_proteins:
                num_edges = self.determine_num_edges_to_cluster(protein, cluster_list, max_edges_until_return=min_num_connections)
                
                if (num_edges >= min_num_connections):
                        qualifying_proteins.append(protein)
                        if (len(qualifying_proteins) >= max_list_length):
                            return qualifying_proteins
            
        
        return qualifying_proteins
        

    def which_proteins_in_cluster_is_a_protein_connected_to() -> list: # or could be of type numpy.ndarray like 'labels' is
        pass



