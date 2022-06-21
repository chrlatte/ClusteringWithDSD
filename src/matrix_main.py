"""
Author: Charlotte Versavel
Date:   June 2022

                             matrix_main.py

Purpose: a main to test the matrix_class

"""

from matrix_class import *

def main():

    matrix = ProteinMatrix("../data/testing_data/tiny_dream3.txt")

    # print(matrix.get_interaction("PRKCA", "GPSM2"))
    print(matrix)

    list = ["PRKCA", "GPSM2", "SRF", "APROTEIN", "OTHERPROTEIN"]
    submatrix = SubMatrix(list, matrix)

    print(f"{submatrix.get_matrix()}")
    # degree = matrix.find_degree("PRKCA")
    # print(degree)

    # degree_list =  DegreeList(matrix)
    # print(f"Degree list: {degree_list.get_degree_list()}")
    # print(f"{degree_list.get_protein_at_index(0, degree=True)}")
     
    


if __name__ == "__main__":
    main()