#!/cluster/tufts/cowenlab/.envs/denoise/bin/python
import sys
import os
sys.path.append(os.getcwd() + "/src")
from dsd_utils import *
import argparse
import pandas as pd
import numpy as np
from clustering import recursive_clustering
from dsd_utils import compute_X_normalized
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics.pairwise import rbf_kernel
import json

def get_args():
<<<<<<< HEAD
	"""
	python clustering_main.py --network "../data/networks/DREAM_files/dream_3.txt"
	"""
	parser = argparse.ArgumentParser()
	parser.add_argument("--network", help = "Network file")
	parser.add_argument("--output", help = "Output file")
	parser.add_argument("--min-c", help = "Min cluster file", default = 50, type = int)
	parser.add_argument("--max-c", help = "Max cluster file", default = 500, type = int)
	return parser.parse_args()
	
=======
    """
    python clustering_main.py --network "../data/networks/DREAM_files/dream_3.txt"
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--network", "-n", help="Network file", required=True)


    return parser.parse_args()


def main(arg):

    network = arg.network
    df = pd.read_csv(network, sep=" ", header=None)
    nodes = set(list(df[0]) + list(df[1]))
    n = len(nodes)

    # Construct a dictionary where each string is associated with an integer value (nodemap)
    # {key: value}, key datatype = string, value = integer
    # KEY: PRKCA, VAL: 0, KEY: GPSM2, VAL:1, A[0,1] = 1
    A = np.zeros((n, n))
    A[4,5]
    print(df)
    print(f"done")

>>>>>>> components

def main(args):
	network = args.network
	output  = args.output
	if (os.path.exists(f"{output}-rbf_0.5.npy") and
	    os.path.exists(f"{output}.json")):
		A1 = np.load(f"{output}-rbf_0.5.npy")
		with open(f"{output}.json", "r") as of:
			nodemap = json.load(of)

	else:	
		df = pd.read_csv(network, sep = " ", header = None)
		nodes = set(list(df[0]) + list(df[1]))

		# get nodemap
		nodemap = {k: i for i, k in enumerate(nodes)}
		with open(f"{output}.json", "w") as of:
			json.dump(nodemap, of)

		print("Constructing the Adjacency matrix...")
		n = len(nodes)
		A = np.zeros((n, n))
		for p, q, w in df.values:
			A[nodemap[p], nodemap[q]] = 1 # make this unweighted
			A[nodemap[q], nodemap[p]] = 1
		
		print("Start Clustering Sequence....")
	    
		# Get the DSD embedding
		print("\tComputing the embedding...")
		if (os.path.exists(f"{output}-dse.npy")):
			X = np.load(f"{output}-dse.npy")
		else:
			X = compute_X_normalized(A)
			np.save(f"{output}-dse.npy", X)

		print("\tComputing the distance and rbf affinity matrix...")
		# squareform pdist
		if (os.path.exists(f"{output}-dist.npy")):
			Dst = np.load(f"{output}-dist.npy")
		else:
			Dst = squareform(pdist(X))
			np.save(f"{output}-dist.npy", Dst)
		# rbf kernel
		A1   = rbf_kernel(Dst, gamma = 0.5)

	print("\t[+]Constructed the new similarity matrix using DSD.")
	print("\tStarting clustering...")
	clusters = recursive_clustering(A1, 3, args.min_c, args.max_c)

	# Reverse nodemap
	r_nodemap = {i:k for k, i in nodemap.items()}

	clustermap = {r_nodemap[c]: label for c, label in enumerate(cluster)}

	with open(f"{output}-cluster.json", "w") as of:
		json.dump(clustermap, of)
	print(clusters)
    
if __name__ == "__main__":
	main(get_args())
