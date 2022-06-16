import pandas as pd 
import numpy as np

class ProteinMatrix:



    # name = "pizza"

	# type = "snack"

    # def _init_(self, type):

   	# 	self.type = type
    
    """
    member variables:
    """

    protein_matrix_df = pd.DataFrame
    proteins_in_matrix = np.array
    protein_indexes = dict()



    def __init__(self, csv_filename):
        # self.type = #dataframe?
        # read in genes and interactions into a dataframe
        
        # print(f"filename: {csv_filename}")

        self.protein_matrix_df = pd.read_csv(csv_filename, delimiter = '\s+', names = ["gene_1", "gene_2", "edge_weight"]) 

        self.proteins_in_matrix = np.unique(np.append(self.protein_matrix_df["gene_1"], self.protein_matrix_df["gene_2"]))

        self._init_dict_of_proteins_and_indexes()
        
        #print(self.protein_matrix_df)
    
    

    def __repr__(self, matrix=False, dict=False): 
        if (not matrix) and (not dict) :
            return "specify <matrix=True> and/or <dict=true> to print an instance of the ProteinMatrix class"
        
        if (matrix):
            with pd.option_context('display.max_rows', 10,
                        'display.max_columns', 10,
                        'display.precision', 5):
                print(self.protein_matrix_df)
        if (dict):
            print(self.protein_indexes)
        
        return ""




    """             _init_dict_of_proteins_and_indexes(self)
        Parameters: none
        Purpose:    a helper function to populate a the protein's index 
                    dictionary 
        Returns:    n/a
    """
    def _init_dict_of_proteins_and_indexes(self):

        index = 0
        for protein in self.proteins_in_matrix:
            self.protein_indexes[protein] = index
            index += 1
        
        