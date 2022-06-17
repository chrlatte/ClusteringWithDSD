import pandas as pd 
import numpy as np

class ProteinMatrix:

    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    * * * * * * * * * * * * * MEMBER VARIABLES * * * * * * * * * * * * * *  
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    protein_data_df = pd.DataFrame
    list_of_all_proteins_in_matrix = np.array
    protein_indexes = dict()
    protein_matrix = pd.DataFrame()


    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    * * * * * * * * * * * * * * INITIALIZERS * * * * * * * * * * * * * * *  
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    def __init__(self, csv_filename: str, **kwargs):

        self.protein_data_df = pd.read_csv(csv_filename, delimiter = '\s+', names = ["gene_1", "gene_2", "edge_weight"]) 

        self.list_of_all_proteins_in_matrix = np.unique(np.append(self.protein_data_df["gene_1"], self.protein_data_df["gene_2"]))
        self._init_dict_of_proteins_and_indexes()

        self._init_matrix()


    def __repr__(self): 
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
            Purpose:    a helper function to populate the matrix with 
                        interactions from a csv file
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



    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    * * * * * * * * * * * * * * * GETTERS * * * * * * * * * * * * * * * * *  
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    def get_matrix(self) -> pd.DataFrame:
        return self.protein_data_df
    
    def get_index(self,  protein : str or None = None) -> dict() or int:
        if (protein == None):
            return self.protein_indexes

        return self.protein_indexes[protein]
    
    def get_protein(self, index : int or None = None) -> np.array or str:
        if (index == None):
            return self.list_of_all_proteins_in_matrix
        
        return self.list_of_all_proteins_in_matrix[index]
        
    def get_interaction(self, index1: int, index2: int):
        return self.protein_matrix.iloc[index1, index2]

class SubMatrix:

    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    * * * * * * * * * * * * * MEMBER VARIABLES * * * * * * * * * * * * * *  
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    list_of_all_proteins_in_matrix = np.array
    protein_indexes = dict() # { 'protein name' : [new_index, old_index]}
    protein_matrix = pd.DataFrame()

    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    * * * * * * * * * * * * * * INITIALIZERS * * * * * * * * * * * * * * *  
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    def __init__(self, proteins_to_map : np.array, original_matrix : ProteinMatrix):
        print(f"instance of submatrix class")
        # note: not all proteins in list of proteins could be in ProteinMatrix. if so, (their index isnt found in the matrix dictionary),catch the error and put the protein in the table but with only zeros. 

        original_index_dict = original_matrix.get_index()
        # original_protein_list = original_matrix.get_protein()

        self.list_of_all_proteins_in_matrix = np.unique(proteins_to_map)
        
        print(self.list_of_all_proteins_in_matrix)


        self._init_dict_of_proteins_and_indexes(original_index_dict)

        self._init_matrix(original_matrix)


        

    

    def _init_dict_of_proteins_and_indexes(self, original_index_dict):
        """             
        Purpose:    a helper function to populate a the protein's index 
                    dictionary 
        Returns:    n/a
        """
        new_index = 0
        for protein in self.list_of_all_proteins_in_matrix:
            try: 
                og_index = original_index_dict[protein]
            except KeyError:
                og_index = -1

            # { "protein name" : [new_index, old_index]}
            self.protein_indexes[protein] = [new_index, og_index]
            
            new_index += 1
        
        print(f"new dict init: {self.protein_indexes}")


    def _init_matrix(self, original_matrix: ProteinMatrix):
        """             
            Purpose:    a helper function to populate the new matrix with 
                        interactions from the original matrix
            Returns:    n/a
        """
        self.protein_matrix = pd.DataFrame(
            columns=self.list_of_all_proteins_in_matrix, 
            index=self.list_of_all_proteins_in_matrix)
        
        self.protein_matrix.fillna(0, inplace=True)
        


        for i in range(np.size(self.list_of_all_proteins_in_matrix)):
            
            protein1 : str = self.list_of_all_proteins_in_matrix[i]

            if (self.protein_indexes[protein1][1] != -1):

                for j in range (i + 1, np.size(self.list_of_all_proteins_in_matrix)) :

                    protein2 : str = self.list_of_all_proteins_in_matrix[j]
                    if (self.protein_indexes[protein2][1] != -1):
                        
                        
                        print(f"interaction between {protein2} and {protein2}: in og matrix")
                    



        # # take each pair of proteins, find their indexes, and then populate the matrix with their interaction
        # for n in range(len(self.protein_data_df)): 
        #     index1 = self.protein_indexes[self.protein_data_df.iloc[n, 0]]
        #     index2 = self.protein_indexes[self.protein_data_df.iloc[n, 1]]
        #     num = self.protein_data_df.iloc[n, 2]

        #     self.protein_matrix.iloc[index1, index2] = num
        #     self.protein_matrix.iloc[index2, index1] = num




# a : int | str = "hello"

# class Foo:
#     def __init__(self, a : int | None = None, b : str | None = None, c : list[str] | None = None) -> None:
#         pass

# f = Foo(b="hello")

# def bar(b, *args, a=3, **kwargs):
#     kwargs["whatever"]
#     ...

# bar(1, 2,3,4, a=18, whatever=3)