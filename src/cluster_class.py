"""
Author: Charlotte Versavel
Date:   June 2022

                             cluster_class.py

Purpose: a class to represent a cluster of proteins

TODO: code documentation / fxn contracts

"""

import pandas as pd 
import numpy as np


class SingleCluster:
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    * * * * * * * * * * * * * MEMBER VARIABLES * * * * * * * * * * * * * *  
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    cluster_number : int = -1
    other_number : float = -1.0
    proteins_in_cluster : np.ndarray = np.array([], dtype=str)

    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    * * * * * * * * * * * * * * INITIALIZERS * * * * * * * * * * * * * * *  
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    def __init__(self, proteins_in_cluster : np.ndarray):
        """             
        Parameters: proteins_in_cluster is an np.ndarray where element 0 is the 
                    cluster number, elem 1 is a float, and the rest of the 
                    array contains protein names 
        Purpose:    to initialize a cluster object from a row in the list of 
                    clusters
        Returns:    n/a
        
        TODO: if the first item is a string, then it is a protein name. if it is an int, then it is a cluster number. may only want to initialize based on whats in the list
        """
        try: 
            self.cluster_number = proteins_in_cluster.pop(0)
            self.other_number = proteins_in_cluster.pop(0)

            self.proteins_in_cluster = proteins_in_cluster
        except:
            print(f"Error while initializing a cluster from array: {proteins_in_cluster}")
        
    def __repr__(self): 
        """             
        Purpose:    to override the print function for this class
        Returns:    a string representation of the cluster
        """
        return f"Cluster {self.cluster_number} has {len(self.proteins_in_cluster)} proteins: {self.proteins_in_cluster}"
        

        
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    * * * * * * * * * * * * * * * GETTERS * * * * * * * * * * * * * * * * *  
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    def get_cluster_number(self) -> int:
        """             
        Purpose:    to allow access to the cluster number
        Returns:    the cluster number
        """
        return int(self.cluster_number)
    
    def get_protein_list(self) -> np.ndarray:
        """             
        Purpose:    to allow access to the list of proteins in a cluster
        Returns:    the list of proteins in the cluster
        """
        return self.proteins_in_cluster
    


class ProteinClusters:

    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    * * * * * * * * * * * * * MEMBER VARIABLES * * * * * * * * * * * * * *  
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    # protein_data_df = pd.DataFrame
    # list_of_all_proteins_in_matrix = np.array
    # protein_indexes = dict()
    # protein_matrix = pd.DataFrame()
    
    all_clusters = np.array([], dtype=SingleCluster)


    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    * * * * * * * * * * * * * * INITIALIZERS * * * * * * * * * * * * * * *  
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    def __init__(self, csv_filename: str, **kwargs):
        """  
        Parameters: csv_filename is the name of a csv file containing several 
                    clusters of proteins           
        Purpose:    to populate several single clusters with data from a CSV 
                    file
        Returns:    n/a
        """
        try:
            # for each line in the file, get that line, and split it into a 
            # list. give that list to the SingleCluster constructor. add that 
            # singlecluster to the list of all clusters
            with open(csv_filename, "r") as data:
                for line in data:
                    line = line.strip().split("\t")
                    # cluster = SingleCluster(line)
                    self.all_clusters = np.append(self.all_clusters, SingleCluster(line))

        except FileNotFoundError:
            print(f"ERROR! file: {csv_filename} not found. Instance of ProteinClusters could not be initialized")

    def __repr__(self): 
        """             
        Purpose:    TODO
        Returns:    TODO
        """
        return f"ProteinClusters has {len(self.all_clusters)} clusters (use the print_all method to see them)"

    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    * * * * * * * * * * * * * * * GETTERS * * * * * * * * * * * * * * * * *  
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    def get_cluster(self, cluster_number: int) -> SingleCluster:
        """             
        Parameters: cluster_number is the number of the cluster to get
        Purpose:    to get a cluster from the list of all clusters
        Returns:    the cluster with the given number
        """
        
        # print(f"get cluster called. printing self.all_clusters: {self.all_clusters}")

        # print(f"looking for cluster number: {cluster_number}")
        # print(f"self.all_clusters[cluster_number] is {self.all_clusters[cluster_number]}. it's cluster number is {(self.all_clusters[cluster_number]).get_cluster_number()}")


        #print(self.all_clusters[cluster_number].get_cluster_number())
        # print(f"comparing {type(cluster_number)} to {type(self.all_clusters[cluster_number].get_cluster_number())}")
        try:
            if (self.all_clusters[cluster_number].get_cluster_number() == cluster_number):
            # first: assume clusters were read in in order
                return self.all_clusters[cluster_number]
            else:
            # possible that clusters were read out of order:
                for cluster in self.all_clusters:
                    if (cluster.get_cluster_number() == cluster_number):
                        return cluster
        except IndexError:
            print(f"ERROR! cluster number {cluster_number} out of bounds. searching for it anyways")
            for cluster in self.all_clusters:
                if (cluster.get_cluster_number() == cluster_number):
                    return cluster
    
        
        print(f"ERROR! cluster {cluster_number} not found. returning empty cluster")
        return SingleCluster([]) # if no cluster was found, return an empty cluster
    
    def get_proteins_from_cluster(self, cluster_number: int) -> np.ndarray:
        """             
        Parameters: cluster_number is the number of the cluster to get
        Purpose:    to get the list of proteins from a cluster
        Returns:    the list of proteins in the cluster
        """
        return self.get_cluster(cluster_number).get_protein_list()

    def get_num_clusters(self) -> int:
        """
        Purpose:    to determine the number of clusters
        Returns:    
        """
        return len(self.all_clusters)


    def print_all(self) -> None:
        """             
        Purpose:    to print all the clusters in the ProteinClusters object
        Returns:    n/a
        """
        for cluster in self.all_clusters:
            print(cluster)