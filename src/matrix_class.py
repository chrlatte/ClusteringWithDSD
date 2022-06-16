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

    protein_data_df = pd.DataFrame
    list_of_all_proteins_in_matrix = np.array
    protein_indexes = dict()
    protein_matrix = pd.DataFrame()




    def __init__(self, csv_filename: str):

        self.protein_data_df = pd.read_csv(csv_filename, delimiter = '\s+', names = ["gene_1", "gene_2", "edge_weight"]) 

        self.list_of_all_proteins_in_matrix = np.unique(np.append(self.protein_data_df["gene_1"], self.protein_data_df["gene_2"]))
        self._init_dict_of_proteins_and_indexes()

        self._init_matrix()


        
    
    

    def __repr__(self, **kwargs): 
        with pd.option_context('display.max_rows', 10,
                        'display.max_columns', 10,
                        'display.precision', 5):
                print(self.protein_matrix)
                # print(self.protein_data_df)
        return ""





    def _init_dict_of_proteins_and_indexes(self):
        """             
        Purpose:    a helper function to populate a the protein's index 
                    dictionary 
        Returns:    n/a
        """
        index = 0
        for protein in self.list_of_all_proteins_in_matrix:
            self.protein_indexes[protein] = index
            index += 1
        



    def _init_matrix(self):
        """             
            Purpose:    a helper function to populate the matrix with interactions from a csv file
            Returns:    n/a
        """
        self.protein_matrix = pd.DataFrame(
            columns=self.list_of_all_proteins_in_matrix, 
            index=self.list_of_all_proteins_in_matrix)
        
        self.protein_matrix.fillna(0, inplace=True)
        
        # take each pair of proteins, find their indexes, and then populate the matrix with their interaction
        for n in range(len(self.protein_data_df)): 
            index1 = self.protein_indexes[self.protein_data_df.iloc[n, 0]]
            index2 = self.protein_indexes[self.protein_data_df.iloc[n, 1]]
            num = self.protein_data_df.iloc[n, 2]

            self.protein_matrix.iloc[index1, index2] = num
            self.protein_matrix.iloc[index2, index1] = num

        

            
            
