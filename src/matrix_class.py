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

    protein_indexes = dict()
    protein_matrix_df = pd.DataFrame


    def __init__(self, csv_filename):
        # self.type = #dataframe?
        # read in genes and interactions into a dataframe
        
        # print(f"filename: {csv_filename}")

        self.protein_matrix_df = pd.read_csv(csv_filename, delimiter = '\s+', names = ["gene_1", "gene_2", "edge_weight"]) 

        #print(self.protein_matrix_df)
    
    

    def __str__(self): 
        return "From str method of Test: a is % s, " \ 
              "b is % s" % (self.a, self.b) 

    def print(*args, **kwargs):
        __builtin__.print('New print function')
        return __builtin__.print(*args, **kwargs)
    