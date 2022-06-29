
"""
Author: Charlotte Versavel
Date:   June 2022

                            cytoscape_print_class.py

Purpose: a class TODO

"""

import pandas as pd 
import numpy as np

from matrix_class import *

class PrintToFile:
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    * * * * * * * * * * * * * MEMBER VARIABLES * * * * * * * * * * * * * *  
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    file: str
    curr_color: str


    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    * * * * * * * * * * * * * * INITIALIZERS * * * * * * * * * * * * * * *  
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    def __init__(self, filepath: str = "cytoscape_input.sif") -> None:
        """            
        Parameters: 
        Purpose:    
        Returns:    
        """
        self.file = filepath
        self.curr_color = 0

        

    def clear_file(self) -> None:
        """
        will all text from the file
        """
        with open(self.file, 'w'):
            pass

    

    def swap_file(self, new_filepath: str) -> None:
        """
        will swap the filepath so future text is printed to the new file
        """
        self.file = new_filepath


    
    def print_matrix(self, matrix: ProteinMatrix or SubMatrix) -> None:
        """
        takes in a matrix, and prints interactions to file in form protein1   protein2, so that the file can be uploaded 
        """
        with open(self.file, 'a') as f:

            for protein in matrix.get_list_of_proteins():
                f.write(f"{protein}")
                connected_proteins = matrix.get_list_of_proteins_connected_to_protein(protein)
                if (connected_proteins): 
                    f.write(f" pp")
                    for connected_protein in connected_proteins:   
                        f.write(f" {connected_protein}")
                f.write("\n")
        

    def assign_colors_to_clusters(self, proteins_in_cluster: list()) -> None:
        """
        takes in a matrix, and prints interactions to file
        """
        with open(self.file, 'a') as f:
            f.write(f"protein: color\n")
            for protein in proteins_in_cluster:
                f.write(f"{protein} {self.curr_color}\n")
            
        
        self.curr_color += 1





        
    