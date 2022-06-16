"""
Author: Charlotte Versavel
Date:   June 2022

                             matrix_main.py

Purpose: a main to test the matrix_class

"""

from matrix_class import *

def main():

    matrix = ProteinMatrix("../data/testing_data/tiny_dream3.txt")

    print(matrix)




if __name__ == "__main__":
    main()