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

        #print(self.protein_matrix_df)
    
    

    def __repr__(self): 
        # return(self.protein_matrix_df.head(10).to_string())

        with pd.option_context('display.max_rows', 10,
                       'display.max_columns', 10,
                       'display.precision', 5):
            print(self.protein_matrix_df)
            return "" # todo: self.protein_matrix_df.to_string() does not 
                      #       return the simplified table
            