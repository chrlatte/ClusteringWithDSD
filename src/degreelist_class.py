"""
Author: Charlotte Versavel
Date:   June 2022

                             degreelist_class.py

Purpose: a class to 

"""


import pandas as pd 
import numpy as np



class DegreeList:
    
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    * * * * * * * * * * * * * MEMBER VARIABLES * * * * * * * * * * * * * *  
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    sorted_protein_degree_dict = dict()
    "* * * * * * * * * * degree dictionary variables * * * * * * * * * * *"
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
        # self.all_proteins = matrix.get_list_of_proteins()

        protein_degree_dict = {name:matrix.find_degree(name) for name in matrix.get_list_of_proteins()}
        # print(f"{protein_degree_dict}")

        self.sorted_protein_degree_dict = sorted(protein_degree_dict.items(), key=lambda x: x[1], reverse=True)

        #print(f"{self.sorted_protein_degree_dict}")

        #print(self.sorted_protein_degree_dict[1])

    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    * * * * * * * * * * * * * * * GETTERS * * * * * * * * * * * * * * * * *  
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    def get_degree_list(self) -> list():
        """             
        Purpose:    to allow access to the sorted degree list
        Returns:    a list of tuples of (protein, degree)
        """
        return self.sorted_protein_degree_dict
    
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
    

    "* * * * * * * * * * degree dictionary functions * * * * * * * * * * *"
    def _init_degree_list_(self) -> None:
        """            
        Purpose:    to take in a proteinMatrix (or submatrix) and create a 
                    sorted dictionary of protein:degree for all proteins in the 
                    matrix.
        Returns:    n/a
        """
        protein_degree_dict = {name:self.find_degree(name) for name in self.list_of_all_proteins_in_matrix}

        self.sorted_protein_degree_dict = sorted(protein_degree_dict.items(), key=lambda x: x[1], reverse=True)


    def get_degree_list(self) -> list():
        """             
        Purpose:    to allow access to the sorted degree list
        Returns:    a list of tuples of (protein, degree)
        """
        return self.sorted_protein_degree_dict
    
    def get_protein_at_degreelist_index(self, index : int, degree = False) -> str or tuple:
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
    
    def create_list_of_proteins_connected_to_cluster(self, cluster_list : np.ndarray, max_list_length : int or None = None, min_num_connections : int = 3) -> np.ndarray:
        """             
        Parameters: cluster_list is a list of proteins in a cluster
                    max_list_length is an upper bound for the number of proteins to return in a list. If None, all proteins with at least min_num_connections connections are added to the list
                    min_num_connections is the minimum number of connections a protein must have to be added to the list and considered 'connected' to the cluster
        Purpose:    to create a list of proteins that are connected to the cluster
        Returns:    a list of proteins that are connected to the cluster
        """
        pass

