"""
Author: Charlotte Versavel
Date:   June 2022

                       json_to_cluster_conversion_main.py

Purpose:    to prototype functions that will take a json file thats a   
            dictionary of proteins and their cluster numbers, and will convert 
            that to a file similar to the other cluster files. 

"""

from matrix_class import *
from cluster_class import *
from degreelist_class import *

import json 


def main():

    clusters = AllClusters()

    with open("../data/testing_data/cluster_dict.txt","r") as cluster_dict_file:#r - open file in read mode
        pass
        dict_of_clusters = json.load(cluster_dict_file)
        # print(dict_of_clusters)
        clusters = AllClusters(protein_to_cluster_dict=dict_of_clusters)

        


    clusters.print_all()




if __name__ == "__main__":
    main()
