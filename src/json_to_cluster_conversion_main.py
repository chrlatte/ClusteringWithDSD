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

    clusters = AllClusters("")

    with open("../data/testing_data/cluster_dict.txt","r") as cluster_dict_file:#r - open file in read mode
        dict_of_clusters = json.load(cluster_dict_file)
        print(f"type: {type(dict_of_clusters)} data: {dict_of_clusters}") #data is a dict

        key_list = list([key for key in dict_of_clusters])

        #print(key_list)

        for protein in key_list:
            clusters.add_protein_to_cluster(protein, dict_of_clusters[protein])


    clusters.print_all()



    with open("../data/results/DREAM-3-cc/d3_5_100_cluster.txt", "w") as cluster_txt_file:

        cluster_txt_file.write('Hi')
        
    
    # # initializing string  
    # string1 = '{"Kiprono" : 67, "Bob" : 76, "Alice" : 88}' 
    
    # # printing original string  
    # print(string1)
    # print ("Type before converting:",type(string)) #check type
    
    # # convert dictionary string to dictionary 
    # res = json.loads(string1)
    
    # # print result 
    # print(res) 
    # print("Type after converting: ",type(res)) #check type
    # ----------Output----------
    # {"Kiprono" : 67, "Bob" : 76, "Alice" : 88}
    # Type before converting: <class 'str'>
    # {'Kiprono': 67, 'Bob': 76, 'Alice': 88}
    # Type after converting:  <class 'dict'>
    # pass


if __name__ == "__main__":
    main()
