"""
Author: Charlotte Versavel
Date:   July 2022

                            analysis_utils.py

Purpose: TODO

"""

from matrix_class import ProteinMatrix
from matrix_class import SubMatrix
from cluster_class import AllClusters
from degreelist_class import DegreeList
import json
from json import load as jsonload


import numpy as np

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
 * * * * * * * * * * * * * * * FUNCTIONS * * * * * * * * * * * * * * *
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
def initialize_matrix_clusters_degreelist(interactions_filepath: str, clusters_filepath: str):
    """
    TODO
    a file that has interactions of the form protein1 TAB protein2 TAB interaction
    a file that containing a dictionary, for each protein in a cluster, the protein's name is linked to its cluster number
    """
    clusters_dict = {}
    # convert actual cluster file to a dictionary!!
    with open(clusters_filepath,"r") as cluster_dict_file:
        clusters_dict = json.load(cluster_dict_file)
    
    matrix = ProteinMatrix(interactions_filepath)
    clusters = AllClusters(protein_to_cluster_dict=clusters_dict)
    degreelist = DegreeList(matrix)

    return matrix, clusters, degreelist



def create_term_mapping_list(go_terms_filepath: str, term_mapping_filepath: str = 'term_mapping.txt'):
    """
    the original file (go_terms_filepath) is in form GOTERM tab PROTEIN, while the term mapping file (term_mapping_filepath) is printed in form PROTEIN tab GOTERM. if a protein has multiple they appear on seperate lines
    """
    with open(term_mapping_filepath, 'w') as file:
        with open(go_terms_filepath, 'r') as go_annotation_file:
            for _ in range(1): # first line of file has column titles, and should be skipped
                next(go_annotation_file)
            for line in go_annotation_file:
                terms = line.split()
                file.write(f"{terms[1]}\t{terms[0]}\n")