#!/usr/bin/python3
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
	"""
	python clustering_main.py --network "../data/networks/DREAM_files/dream_3.txt"
	"""
	parser = argparse.ArgumentParser()
	parser.add_argument("--network", help = "Network file")
	parser.add_argument("--output", help = "Output file")
	parser.add_argument("--min-c", help = "Min cluster file", default = 3, type = int)
	parser.add_argument("--max-c", help = "Max cluster file", default = 500, type = int)
	parser.add_argument("--spectral", help = "Spectral clustering parameter", default = 50, type = int)
	parser.add_argument("--simple-inverse", help = "Use simple inverse instead of RBF kernal")
	parser.add_argument("--remove-low-weights", help = "Remove low weight scores in the RBF similarity matrix")
	parser.add_argument("--seed", help = "Random seed for spectral clustering parameter", default = 68, type = int)
	return parser.parse_args()
	

def main(args):
	print(f"Starting clustering with the following parameters")
	print(f"\tNetwork: {args.network}")
	print(f"\tMinimim cluster size: {args.min_c}")
	print(f"\tMaximim cluster size: {args.max_c}")
	print(f"\tSpectral parameter: {args.spectral}")
	sim_metric = "Simple inverse" if args.simple_inverse else "RBF"
	print(f"\tSimilarity metric: {sim_metric}")
	network = args.network
	output  = args.output
	sim_fp  = f"{output}-si.npy" if args.simple_inverse else f"{output}-rbf_0.5.npy"
	if (os.path.exists(f"{output}-rbf_0.5.npy") and
	    os.path.exists(f"{output}.json")):
		A1 = np.load(sim_fp)
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
		if args.simple_inverse:
			# use simple inverse
			sim_fn = np.vectorize(lambda x: 1 if x == 0 else 1 / float(x))
			A1 = sim_fn(Dst)
			np.save(sim_fp, A1)
		else:
			# rbf kernel
			A1   = rbf_kernel(Dst, gamma = 0.5)
			if args.remove_low_weights:
				A1   = np.where(A1 > 0.0005, A1, 0)
			np.save(sim_fp, A1)

	print("\t[+]Constructed the new similarity matrix using DSD.")
	print("\tStarting clustering...")
	clusters = recursive_clustering(A1, args.spectral, args.min_c, args.max_c, args.seed)

	# Reverse nodemap
	r_nodemap = {i:k for k, i in nodemap.items()}

	clustermap = {}
	for cluster, labels in clusters.items():
		for label in labels:
			clustermap[r_nodemap[label]] = cluster
	with open(f"{output}-cluster1.json", "w") as of:
		json.dump(clustermap, of)

	clustermap = {r_nodemap[c]: label for c, label in enumerate(clusters)}

	with open(f"{output}-cluster.json", "w") as of:
		json.dump(clustermap, of)
	print(clusters)
    
if __name__ == "__main__":
	main(get_args())
