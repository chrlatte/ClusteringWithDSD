
"""
Author: Charlotte Versavel
Date:   June 2022

                            cytoscape_print_class.py

Purpose: a class TODO

"""

from cluster_class import AllClusters

from matrix_class import *

class PrintToFile:
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    * * * * * * * * * * * * * MEMBER VARIABLES * * * * * * * * * * * * * *  
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    * * * * * * * * * * * * * * INITIALIZERS * * * * * * * * * * * * * * *  
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    def __init__(self) -> None:
        """            
        Parameters: 
        Purpose:    
        Returns:    
        """
        pass
    
    
    def print_all_interactions(self, matrix: ProteinMatrix or SubMatrix, filepath:str = "cytoscape_input.sif") -> None:
        """
        Parameters:
            -   filepath is the path of the file to print the interactions to. 
                it defaults to "cytoscape_input.sif" but can be changed to a 
                specified value
        Purpose:    takes in a matrix, and prints interactions to file in a SIF 
                    form (protein1 <interaction> protein2 protein3 ...) so that the file can be uploaded as a network to cytoscape
        Returns: n/a
        """
        with open(filepath, 'w') as f:

            for protein in matrix.get_list_of_proteins():
                f.write(f"{protein}")

                connected_proteins = matrix.get_list_of_proteins_connected_to_protein(protein)

                if (connected_proteins): 
                    f.write(f" pp")
                    for connected_protein in connected_proteins:   
                        f.write(f" {connected_protein}")
                
                f.write("\n")
        

    def assign_colors_to_clusters(self, clusters: AllClusters, filepath:str = "cytoscape_colors.txt") -> None:
        """
        Parameters:
            -   filepath is the path of the file to print the colors to. it 
                defaults to "cytoscape_colors.txt" but can be changed
        TODO
        """
        with open(filepath, 'w') as f:
            f.write(f"protein:\tcolor\n")

            for i in range(clusters.get_num_clusters()):
                
                proteins = clusters.get_cluster_proteins(i)
                for protein in proteins:
                    f.write(f"{protein}\t{i}\n")
            
    
    def print_all_proteins(self, matrix: ProteinMatrix or SubMatrix, filepath:str = "protein_list.txt") -> None:
        
        with open(filepath, 'w') as f:
            
            for protein in matrix.get_list_of_proteins():
                f.write(f"{protein}\n")





        
    