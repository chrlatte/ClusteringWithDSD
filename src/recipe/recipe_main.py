"""
Author: Charlotte Versavel
Date:   July 2022

                            recipe_main.py

Purpose: a main TODO

"""
import numpy as np
import pandas as pd
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import argparse

from classes.matrix_class import ProteinMatrix
from classes.cluster_class import AllClusters
from classes.degreelist_class import DegreeList

from func_e.FUNC_E import FUNC_E 
import func_e.vocabs.all as vocabs


from recipe_utils import initialize_matrix_clusters_degreelist
from recipe_utils import find_clusters_and_proteins_together
# from recipe_utils import pick_ratio
from recipe_utils import print_querylist_of_clusters_to_file
from recipe_utils import get_initialized_fe
from recipe_utils import print_protein_background_to_file
from recipe_utils import *




def get_args():
    """
    TODO
    """
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--network", help = "[required] network file in form Protein1 TAB Protein2 TAB Interaction", default = "data/dream_3.txt", type=str)
    
    parser.add_argument("--cluster_dict", help = "[required] please include a file that is a dictionary of clusters, linking protein name to cluster num ( { protein_name : cluster_num } ) ", default = "data/protein_to_cluster_dict.json")

    parser.add_argument("--protein_to_goterm", help = "[required] please include a file that maps each protein to it's corresponding GO term in form protein TAB GO term", default = "data/term_mapping.txt")

    parser.add_argument("--protein_list", help = "[optional] please include a filepath for a file containing all proteins in the network that you are using. Alternatively, you can choose to not include this, and one will be generated from the given network", default= "None", type=str)

    parser.add_argument("--original_querylist", help = "[optional] if you have a preference for their names, specify a name for the querylist files", default= "data/original_querylist.txt", type=str)
    parser.add_argument("--updated_querylist", help = "[optional] if you have a preference for their names, specify a name for the querylist files", default= "data/updated_querylist.txt", type=str)
    



	# parser.add_argument("--output", help = "Output file")
	# parser.add_argument("--min-c", help = "Min cluster file", default = 50, type = int)
	# parser.add_argument("--max-c", help = "Max cluster file", default = 500, type = int)
    return parser.parse_args()




def main(args):

    print(f"getting command line arguments")

    # setting up files:
    cluster_dict_filepath = args.cluster_dict
    network_filepath = args.network
    term_mapping_filepath = args.protein_to_goterm

    original_querylist_filepath = args.original_querylist
    updated_querylist_filepath = args.updated_querylist

    protein_list_filepath = args.protein_list

    termlist = vocabs.getTerms(['GO'])


    # read data into data structures
    matrix, clusters, degreelist = initialize_matrix_clusters_degreelist(network_filepath, cluster_dict_filepath)
    

    print(f"finding qualifying clusters and proteins")
    ## TODO: make fxn parameters inputs from command line
    # find qualifying clusters and proteins
    qualifying_clusters, qualifying_proteins = find_clusters_and_proteins_together(matrix, clusters, degreelist, cluster_ratio=0, cluster_constant=5, protein_ratio=1, protein_constant=0, use_sqrt=True)


    # TODO: can have an if conditional: if want to print all results use clusters.get_all_cluster_labels or qualifying_clusters, if want only targeted clusters use qualifying_proteins.keys()
    

    print(f"printing querylists to files")
    # original:
    print_querylist_of_clusters_to_file(clusters, qualifying_proteins.keys(),original_querylist_filepath)
    # clusters with added proteins:
    print_querylist_of_clusters_to_file(clusters, qualifying_proteins.keys(), updated_querylist_filepath, qualifying_proteins)


    if protein_list_filepath == "None":
        print_protein_background_to_file(matrix, filename="data/background_proteinlist.txt")
        protein_list_filepath = "data/background_proteinlist.txt"
    
    print(f"doing functional enrichment")
    original_fe = get_initialized_fe(protein_list_filepath, term_mapping_filepath, termlist = termlist)
    original_fe.importFiles({ 'query': original_querylist_filepath })
    original_fe.run(cluster=False)

    updated_fe = get_initialized_fe(protein_list_filepath, term_mapping_filepath, termlist = termlist)
    updated_fe.importFiles({ 'query': updated_querylist_filepath })
    updated_fe.run(cluster=False)


    print(f"assembling data")
    # original:
    original_df = original_fe.enrichment[['Module', 'Term', 'Fishers_pvalue']].copy()
    original_df['Module'] = original_df['Module'] + ' ' + original_df['Term']
    original_df.drop('Term', axis=1, inplace=True)

    # clusters with added proteins:
    updated_df = updated_fe.enrichment[['Module', 'Term', 'Fishers_pvalue']].copy()

    updated_df['Module'] = updated_df['Module'] + ' ' + updated_df['Term']
    updated_df.drop('Term', axis=1, inplace=True)
    updated_df.rename(columns = {'Fishers_pvalue':'Updated_Fishers_pvalue'}, inplace = True)

    # combine into a single df:
    results_df = pd.merge(original_df, updated_df, on=['Module'], how='outer')
    results_df.plot(x="Module", y=["Fishers_pvalue", "Updated_Fishers_pvalue"], kind="bar")

    print(f"number of significant pvals in original clusters: {results_df['Fishers_pvalue'].count()}")
    print(f"number of significant pvals in updated clusters: {results_df['Updated_Fishers_pvalue'].count()}\n")

    print(f"number of clusters that were updated: {len(qualifying_proteins.keys())}")
    print(f"number of clusters that were functionally enriched: {original_fe.enrichment['Module'].nunique()}")
    print(f"number of clusters that were functionally enriched after being updated: {updated_fe.enrichment['Module'].nunique()}")




if __name__ == "__main__":
    main(get_args())