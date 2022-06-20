"""
Author: Charlotte Versavel
Date:   June 2022

                             cluster_main.py

Purpose: a main to test the cluster class

"""


from cluster_class import *

def main():

    # read in the cluster file into a dataframe
    cluster_df = pd.read_csv("../data/testing_data/cluster_test.txt", sep="\t", header=0)
    # create a cluster object
    cluster = Cluster(cluster_df)

    # print the cluster object
    print(cluster)
    



if __name__ == "__main__":
    main()