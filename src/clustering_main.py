from dsd_utils import *
import argparse
import pandas as pd
import numpy as np
from clustering import recursive_clustering
from dsd_utils import compute_X_normalized
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics.pairwise import rbf_kernel

def get_args():
    """
    python clustering_main.py --network "../data/networks/DREAM_files/dream_3.txt"
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--network", help = "Network file")
    return parser.parse_args()


def main(args):
    network = args.network
    df = pd.read_csv(network, sep = " ", header = None)
    nodes = set(list(df[0]) + list(df[1]))

    # get nodemap
    nodemap = {k: i for i, k in enumerate(nodes)}

    n = len(nodes)
    A = np.zeros((n, n))
    for p, q, w in df.values:
        A[nodemap[p], nodemap[q]] = 1 # make this unweighted
        A[nodemap[q], nodemap[p]] = 1
    print("Start Clustering Sequence....")
        
    # Get the DSD embedding
    print("\tComputing the embedding...")
    X = compute_X_normalized(A)

    print("\tComputing the distance and rbf affinity matrix...")
    # squareform pdist
    Dst = squareform(pdist(X))

    # rbf kernel
    A1   = rbf_kernel(Dst, gamma = 0.5)

    print("\t[+]Constructed the new similarity matrix using DSD.")
    print("\tStarting clustering...")
    clusters = recursive_clustering(A1, 3, 500, 100)
    print(clusters)
    
    
if __name__ == "__main__":
    main(get_args())
