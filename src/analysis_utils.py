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
import pandas as pd
import func_e.vocabs.all as vocabs
from func_e.FUNC_E import FUNC_E 


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



# def print_both_querylists_to_files(qualifying_clusters: list, original_clusters: AllClusters, new_clusters: AllClusters, original_query_filepath:str = 'original_querylist.txt', new_query_filepath:str = 'new_querylist.txt') -> None:
#     """TODO"""
#     original_clusters.print_querylist_of_clusters_to_file(qualifying_clusters, query_filepath=original_query_filepath)
#     new_clusters.print_querylist_of_clusters_to_file(qualifying_clusters, query_filepath=new_query_filepath)


def get_initialized_fe(background_filepath: str, terms2features_filepath: str, termlist: pd.DataFrame() = vocabs.getTerms(['GO']), ecut: float = 0.01) -> FUNC_E():
    """TODO"""
    fe = FUNC_E()

    fe.importFiles({
        'background': background_filepath, 
        'terms2features': terms2features_filepath })
    fe.setTerms(termlist)
    fe.setEnrichmentSettings({'ecut': ecut})

    # now all that is left to do is upload the querylist using fe.importFiles({'query': querylist }), and running it, using fe.run(cluster=False)
    return fe


def print_querylist_of_clusters_to_file(clusters: AllClusters, clusters_to_print: list(), query_filepath: str = "querylist.txt", proteins_to_add: dict() = dict()):
    """
    clusters_to_print -> specify a list of which clusters to print
    TODO
    proteins_to_add -> dictionary containing key: clusternum and value: list of proteins to add to that cluster
    """
    output_file = open(query_filepath, 'w')

    for cluster_num in clusters_to_print:
        for protein in clusters.get_cluster_proteins(cluster_num):
            output_file.write(f"{protein}\tcluster_{cluster_num}\n")
    
    if proteins_to_add: # dict not empty -> dict contains 
        for cluster_num in proteins_to_add: 
            for protein in proteins_to_add[cluster_num]: 
                output_file.write(f"{protein}\tcluster_{cluster_num}\n")

    output_file.close()