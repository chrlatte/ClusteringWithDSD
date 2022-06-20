"""
Author: Charlotte Versavel
Date:   June 2022

                             matrix_class.py

Purpose: a class to represent a matrix of protein interactions

"""


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

    "* * * * * * * * * * degree dictionary variables * * * * * * * * * * *"
    sorted_protein_degree_dict = dict()


    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    * * * * * * * * * * * * * * INITIALIZERS * * * * * * * * * * * * * * *  
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    def __init__(self, csv_filename: str, **kwargs):
        """             
        Purpose:    to populate a matrix with data from a CSV file
        Returns:    n/a
        """
        try:
            # read csv file into a dataframe
            self.protein_data_df = pd.read_csv(csv_filename, delimiter = '\s+', 
                names = ["gene_1", "gene_2", "edge_weight"]) 
            # store names of all proteins
            self.list_of_all_proteins_in_matrix = np.unique(np.append(self.
                protein_data_df["gene_1"], self.protein_data_df["gene_2"]))
            # populate the dictionary relating protein names to their indexes
            self._init_dict_of_proteins_and_indexes()
            # populate the matrix with protein interactions
            self._init_matrix()
        
        except FileNotFoundError:
            print(f"ERROR! file: {csv_filename} not found. ProteinMatrix could not be initialized")


    def __repr__(self): 
        """             
        Purpose:    to override the print function for this class so only a 
                    portion of the matrix is shown
        Returns:    an empty string. all printing is done in the function
        """
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
        Purpose:    a helper function to populate the matrix with interactions
                    from a csv file
        Returns:    n/a
        """
        self.protein_matrix = pd.DataFrame(
            columns=self.list_of_all_proteins_in_matrix, 
            index=self.list_of_all_proteins_in_matrix)
        
        self.protein_matrix.fillna(0, inplace=True)
        
        # take each pair of proteins, find their indexes, and then populate the 
        # matrix with their interaction
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
        """             
        Purpose:    allows for access to a 2D matrix of protein interactions
        Returns:    a dataframe of protein interactions
        """
        return self.protein_data_df
    
    def get_list_of_proteins(self) -> np.array:
        """             
        Purpose:    to access the list of all proteins in the matrix
        Returns:    an array of all proteins in the matrix
        """
        return self.list_of_all_proteins_in_matrix

    def get_index(self,  protein : str or None = None) -> dict() or int:
        """             
        Purpose:    to allow access to the protein->index dictionary
        Returns:    either the index for a specific protein or a dictionary of 
                    protein names and their indexes
        """
        if (protein == None):
            return self.protein_indexes

        return self.protein_indexes[protein]
    
    def get_protein(self, index : int or None = None) -> np.array or str:
        """             
        Purpose:    to access the names of proteins in the matrix
        Returns:    either the name of a protein at a specific index, or the 
                    array of all proteins.
        """
        if (index == None):
            return self.list_of_all_proteins_in_matrix
        
        return self.list_of_all_proteins_in_matrix[index]
        
    def get_interaction(self, index1: int, index2: int):
        """             
        Purpose:    to access the interaction values stored in the matrix
        Returns:    the value at the specified indexes
        """
        return self.protein_matrix.iloc[index1, index2]

    def find_degree(self, protein: str) -> int:
        """             
        Parameters: protein the name of a protein
        Purpose:    to find the degree of a protein
        Returns:    the degree of the protein
        """
        count = 0
        for elem in self.protein_matrix.loc[protein]:
            if elem != 0:
                count += 1
        
        return count



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








# TODO: descriptions of the two classes and their functions
# TODO: test the getters for submatrix



class SubMatrix:

    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    * * * * * * * * * * * * * MEMBER VARIABLES * * * * * * * * * * * * * *  
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    list_of_all_proteins_in_matrix = np.array
    protein_indexes = dict() # { 'protein name' : new_index }
    protein_matrix = pd.DataFrame()

    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    * * * * * * * * * * * * * * INITIALIZERS * * * * * * * * * * * * * * *  
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    def __init__(self, proteins_to_map : np.array, original_matrix : ProteinMatrix):
        """             
        Purpose:    to populate the submatrix with data from the original 
                    matrix. 
        Returns:    n/a
        """
        # initialize list of proteins:
        self.list_of_all_proteins_in_matrix = np.unique(proteins_to_map)
        # inititalize new dictionary for the proteins in this submatrix:
        self._init_dict_of_proteins_and_indexes(original_matrix.get_index())
        # inititalize matrix:
        self._init_matrix(original_matrix)
    

    def _init_dict_of_proteins_and_indexes(self, original_index_dict):
        """             
        Purpose:    a helper function to populate a the protein's index 
                    dictionary in form { 'protein_name' = protein_index }
        Returns:    n/a
        """
        new_index = 0
        for protein in self.list_of_all_proteins_in_matrix:
            
            self.protein_indexes[protein] = new_index
            new_index += 1
        

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

            try: # find protein 1 in original matrix:
                index1 = original_matrix.get_index(protein1)
                # fill in all of protein 1's interactions with other proteins in the submatrix:
                for j in range (i + 1, np.size(self.list_of_all_proteins_in_matrix)):
                    protein2 : str = self.list_of_all_proteins_in_matrix[j]

                    try: # find protein 2 in original matrix:
                        index2 = original_matrix.get_index(protein2)
                        # get the interaction between p1 and p2 from original matrix
                        interaction = original_matrix.get_interaction(index1, index2)
                        # populate submatrix
                        self.protein_matrix.iloc[i, j] = interaction
                        self.protein_matrix.iloc[j, i] = interaction

                    except KeyError: # protein2 not in original matrix
                        pass
            except KeyError: # protein1 not in original matrix
                pass 

    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    * * * * * * * * * * * * * * * GETTERS * * * * * * * * * * * * * * * * *  
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    def get_matrix(self) -> pd.DataFrame:
        """             
        Purpose:    allows for access to a 2D matrix of protein interactions
        Returns:    a dataframe of protein interactions
        """
        return self.protein_data_df
    
    def get_list_of_proteins(self) -> np.array:
        """             
        Purpose:    allows for access to the list of proteins in the matrix
        Returns:    a list of proteins
        """
        return self.list_of_all_proteins_in_matrix
    
    def get_index(self,  protein : str or None = None) -> dict() or int:
        """             
        Purpose:    to allow access to the protein->index dictionary
        Returns:    either the index for a specific protein or a dictionary of 
                    protein names and their indexes
        """
        if (protein == None):
            return self.protein_indexes

        return self.protein_indexes[protein]
    
    def get_protein(self, index : int or None = None) -> np.array or str:
        """             
        Purpose:    to access the names of proteins in the matrix
        Returns:    either the name of a protein at a specific index, or the 
                    array of all proteins.
        """
        if (index == None):
            return self.list_of_all_proteins_in_matrix
        
        return self.list_of_all_proteins_in_matrix[index]
        
    def get_interaction(self, index1: int, index2: int):
        """             
        Purpose:    to access the interaction values stored in the matrix
        Returns:    the value at the specified indexes
        """
        return self.protein_matrix.iloc[index1, index2]

    def find_degree(self, protein: str) -> int:
        """             
        Parameters: protein the name of a protein
        Purpose:    to find the degree of a protein
        Returns:    the degree of the protein
        """
        count = 0
        for elem in self.protein_matrix.loc[protein]:
            if elem != 0:
                count += 1
        
        return count

# a : int | str = "hello"

# class Foo:
#     def __init__(self, a : int | None = None, b : str | None = None, c : list[str] | None = None) -> None:
#         pass

# f = Foo(b="hello")

# def bar(b, *args, a=3, **kwargs):
#     kwargs["whatever"]
#     ...

# bar(1, 2,3,4, a=18, whatever=3)



# class DegreeList:
    
#     """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#     * * * * * * * * * * * * * MEMBER VARIABLES * * * * * * * * * * * * * *  
#     """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#     sorted_protein_degree_dict = dict()

#     """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#     * * * * * * * * * * * * * * INITIALIZERS * * * * * * * * * * * * * * *  
#     """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#     def __init__(self, matrix : ProteinMatrix) -> None:
#         """            
#         Parameters: matrix is populated with proteins and their interaction 
#                     weights
#         Purpose:    to take in a proteinMatrix (or submatrix) and create a 
#                     sorted dictionary of protein:degree for all proteins in the 
#                     matrix.
#         Returns:    n/a
#         """
#         # self.all_proteins = matrix.get_list_of_proteins()

#         protein_degree_dict = {name:matrix.find_degree(name) for name in matrix.get_list_of_proteins()}
#         # print(f"{protein_degree_dict}")

#         self.sorted_protein_degree_dict = sorted(protein_degree_dict.items(), key=lambda x: x[1], reverse=True)

#         #print(f"{self.sorted_protein_degree_dict}")

#         #print(self.sorted_protein_degree_dict[1])

#     """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#     * * * * * * * * * * * * * * * GETTERS * * * * * * * * * * * * * * * * *  
#     """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#     def get_degree_list(self) -> list():
#         """             
#         Purpose:    to allow access to the sorted degree list
#         Returns:    a list of tuples of (protein, degree)
#         """
#         return self.sorted_protein_degree_dict
    
#     def get_protein_at_index(self, index : int, degree = False) -> str or tuple:
#         """             
#         Parameters: index is the index of the protein in the sorted list
#                     degree is a boolean that determines if the degree is returned as well
#         Purpose:    to return the protein at the specified index
#         Returns:    the protein at the specified index, or if degree is True, a 
#                     tuple of (protein, degree)
#         """
#         if not degree:
#             return self.sorted_protein_degree_dict[index][0]
#         else:
#             return self.sorted_protein_degree_dict[index]
    

