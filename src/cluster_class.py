"""
Author: Charlotte Versavel
Date:   June 2022

                             cluster_class.py

Purpose: a class to represent a cluster of proteins

"""

import pandas as pd 
import numpy as np

class ProteinCluster:
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    * * * * * * * * * * * * * MEMBER VARIABLES * * * * * * * * * * * * * *  
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    df_row = pd.Series() # display(df.iloc[n])

    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    * * * * * * * * * * * * * * INITIALIZERS * * * * * * * * * * * * * * *  
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    def __init__(self, row = None):
        """             
        Purpose:    to populate the submatrix with data from the original 
                    matrix. 
        Returns:    n/a
        """
        self.df_row = row
        pass


    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    * * * * * * * * * * * * * * * GETTERS * * * * * * * * * * * * * * * * *  
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    def get_index(self,  protein : str or None = None) -> dict() or int:
        """             
        Purpose:    to allow access to the protein->index dictionary
        Returns:    either the index for a specific protein or a dictionary of 
                    protein names and their indexes
        """
        pass
