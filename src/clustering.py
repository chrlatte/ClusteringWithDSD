from sklearn.cluster import SpectralClustering
import numpy as np
import sys
import os
sys.path.append(os.getcwd())

def recursive_clustering(A, clusters, lower, higher):
    indices_to_cluster = [list(range(A.shape[0]))]
    labels             = {}
    label_count        = 0
    while(True):
        
        # Break condition
        if len(indices_to_cluster) == 0:
            break
        
        ids   = indices_to_cluster.pop()
        Aid   = A[np.ix_(ids, ids)] 
        idmap = {i:k for i, k in enumerate(ids)}
        c_ss  = SpectralClustering(n_clusters = clusters, affinity="precomputed").fit(Aid)
        
        # labels - id mapping
        gen_labels = {}
        for i, label in enumerate(c_ss.labels_):
            gen_labels[label] = ([idmap[i]] if label not in gen_labels
                                 else gen_labels[label] + [idmap[i]])
        for l in gen_labels:
            if len(gen_labels[l]) < lower:
                pass
            elif len(gen_labels[l]) > higher:
                indices_to_cluster.append(gen_labels[l].copy())
            else:
                labels[label_count] = [k for k in gen_labels[l]]
                label_count        += 1
    return labels
