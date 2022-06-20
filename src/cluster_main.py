"""
Author: Charlotte Versavel
Date:   June 2022

                             cluster_main.py

Purpose: a main to test the cluster class

"""


from cluster_class import *

def main():

    # read in the cluster file into a dataframe
    cluster_df = pd.read_csv("../data/testing_data/small_cluster.txt", sep="\t", header=0) # TODO: pandas.errors.ParserError: Error tokenizing data. C error: Expected 39 fields in line 31, saw 92 -> need to read line by line instead of reading in the entire file and have it create arrays of the data
    
    print(cluster_df)

    # create a cluster object
    cluster = ProteinCluster(cluster_df)

    # print the cluster object
    print(cluster)
    



if __name__ == "__main__":
    main()