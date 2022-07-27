"""
Author: Charlotte Versavel
Date:   July 2022

                            analysis_main.py

Purpose: a main TODO

"""

import numpy as np
import pandas as pd
import argparse




from matrix_class import ProteinMatrix
from matrix_class import SubMatrix
from cluster_class import AllClusters
from degreelist_class import DegreeList


from connected_components_utils import pick_ratio
from connected_components_utils import find_clusters_and_proteins_together

from analysis_utils import *

from func_e.FUNC_E import FUNC_E # import the class
import func_e.vocabs.all as vocabs



def get_args():
	"""
	python clustering_main.py --network "../data/networks/DREAM_files/dream_3.txt"
	"""
	parser = argparse.ArgumentParser()
	# parser.add_argument("--cluster_dict", help = "please include a file that is a dictionary of clusters, linking protein name to cluster num ( { protein_name : cluster_num } ) ", default = "../data/results/DREAM-3-cc/d3_5_100.json-cluster.json", type=str)
	# parser.add_argument("--network", help = "Network file", default = "../data/networks/DREAM_files/dream_3.txt", type=str)



	# parser.add_argument("--output", help = "Output file")
	# parser.add_argument("--min-c", help = "Min cluster file", default = 50, type = int)
	# parser.add_argument("--max-c", help = "Max cluster file", default = 500, type = int)

	return parser.parse_args()




def main(args):

    cluster_dict_filepath = args.cluster_dict
    network_filepath = args.network
    list_of_proteins_in_cluster_filepath = ""
    termlist = vocabs.getTerms(['GO'])


    matrix, clusters, degreelist = initialize_matrix_clusters_degreelist(network_filepath, cluster_dict_filepath)
    

    genomic_background_filepath = '../data/testing_data/protein_list.txt'
    term_mapping_filepath = 'term_mapping.txt'
    create_term_mapping_list('../data/go-results/dream_3_go_results.tsv', term_mapping_filepath)




    # qualifying_clusters, qualifying_proteins = find_clusters_and_proteins_together(matrix, clusters, degreelist, cluster_ratio=pick_ratio(clusters.get_num_clusters()), protein_ratio=.05, protein_constant=2)



    # og_clusters_querylist_path = 'original_querylist.txt'
    # clusters.print_querylist_of_clusters_to_file(qualifying_clusters, og_clusters_querylist_path)



    # original_fe = FUNC_E()
    # original_fe.importFiles({
    #     'background': genomic_background_filepath, 
    #     'query': og_clusters_querylist_path, 
    #     'terms2features': term_mapping_filepath })
    # original_fe.setTerms(termlist)
    # original_fe.setEnrichmentSettings({'ecut': 0.01})
    # original_fe.run(cluster=False)



    # update_clusters(clusters, qualifying_proteins)


    # updated_clusters_querylist_path = 'new_querylist.txt'
    # clusters.print_querylist_of_clusters_to_file(qualifying_clusters, updated_clusters_querylist_path)

    # updated_fe = FUNC_E()
    # updated_fe.importFiles({
    #     'background': genomic_background_filepath, 
    #     'query': updated_clusters_querylist_path, 
    #     'terms2features': term_mapping_filepath })
    # updated_fe.setTerms(termlist)
    # updated_fe.setEnrichmentSettings({'ecut': 0.01})
    # updated_fe.run(cluster=False)
    
    # # print(f"ORIGINAL CLUSTERS TABLE:\n{original_fe.enrichment.sort_values('Fishers_pvalue')}")
    # # print(f"UPDATED CLUSTERS TABLE:\n{updated_fe.enrichment.sort_values('Fishers_pvalue')}")


    # updated_df = updated_fe.enrichment[['Module', 'Term', 'Fishers_pvalue']].copy()
    # original_df = original_fe.enrichment[['Module', 'Term', 'Fishers_pvalue']].copy()

    # print(f"ORIGINAL CLUSTERS TABLE:\n{original_df}")
    # print(f"UPDATED CLUSTERS TABLE:\n{updated_df}")

    # results_df = pd.merge(original_df, updated_df, on=['Module','Term'], how='outer')

    # print(f"FINAL MERGE:\n{results_df}")
    
    # results_df.to_csv("results.txt", sep='\t')


if __name__ == "__main__":
    main(get_args())