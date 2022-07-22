"""
Author: Charlotte Versavel
Date:   July 2022

                            analysis_main.py

Purpose: a main TODO

"""

import numpy as np
import json


from matrix_class import ProteinMatrix
from matrix_class import SubMatrix
from cluster_class import AllClusters
from degreelist_class import DegreeList


from connected_components_utils import *
from analysis_utils import *

from func_e.FUNC_E import FUNC_E # import the class
import func_e.vocabs.all as vocabs


def main():

    matrix, clusters, degreelist = initialize_matrix_clusters_degreelist("../data/networks/DREAM_files/dream_3.txt", "../data/results/DREAM-3-cc/d3_5_100.json-cluster.json")
    
    qualifying_clusters, qualifying_proteins = find_clusters_and_proteins_together(matrix, clusters, degreelist, cluster_ratio=pick_ratio(clusters.get_num_clusters()), protein_ratio=.05, protein_constant=2)


    genomic_background_filepath = '../data/testing_data/protein_list.txt'
    termlist = vocabs.getTerms(['GO'])
    term_mapping_filepath = 'term_mapping.txt'
    create_term_mapping_list('../data/go-results/dream_3_go_results.tsv', term_mapping_filepath)



    og_clusters_querylist_path = 'original_querylist.txt'
    clusters.print_querylist_of_clusters_to_file(qualifying_clusters, og_clusters_querylist_path)



    original_fe = FUNC_E()
    original_fe.importFiles({
        'background': genomic_background_filepath, 
        'query': og_clusters_querylist_path, 
        'terms2features': term_mapping_filepath })
    original_fe.setTerms(termlist)
    original_fe.setEnrichmentSettings({'ecut': 0.01})
    original_fe.run(cluster=False)



    update_clusters(clusters, qualifying_proteins)


    updated_clusters_querylist_path = 'new_querylist.txt'
    clusters.print_querylist_of_clusters_to_file(qualifying_clusters, updated_clusters_querylist_path)

    updated_fe = FUNC_E()
    updated_fe.importFiles({
        'background': genomic_background_filepath, 
        'query': updated_clusters_querylist_path, 
        'terms2features': term_mapping_filepath })
    updated_fe.setTerms(termlist)
    updated_fe.setEnrichmentSettings({'ecut': 0.01})
    updated_fe.run(cluster=False)
    
    print(f"ORIGINAL CLUSTERS TABLE:\n{original_fe.enrichment.sort_values('Fishers_pvalue')}")
    print(f"UPDATED CLUSTERS TABLE:\n{updated_fe.enrichment.sort_values('Fishers_pvalue')}")
    

if __name__ == "__main__":
    main()
