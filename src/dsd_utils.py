# Requirements: Install numpy and scipy using 
#       pip install numpy scipy

# If you are in a conda environemnt, activate the environment
# and type 
#       conda install numpy scipy

import numpy as np
import scipy as sp
from scipy.spatial.distance import pdist, squareform

def compute_X_normalized(A, 
                        t = -1, 
                        lm = 1, 
                        is_normalized = False):
    """
    Computes the DSD embedding from the adjacency matrix "A"
    In order to use the cDSD, basically set `t=-1`, `lm=1` and `is_normalized=False`
    """
    
    # compute the degree vector
    d   = A @ np.ones((A.shape[0], 1))

    # To ensure a small value at location where d == 0,
    d   = np.where(d > 0, d, 0.00001)

    P   = A / d

    Identity = np.identity(A.shape[0])
    e = np.ones((A.shape[0], 1))

    # Compute W
    scale = np.matmul(e.T, d * e)[0, 0]
    W = np.multiply(1 / scale, np.matmul(e, (d * e).T))

    up_P = np.multiply(lm, P - W)
    X_ = Identity - up_P
    X_i = np.linalg.pinv(X_)

    if t > 0:
        LP_t = Identity - np.linalg.matrix_power(up_P, t)
        X_i = np.matmul(X_i, LP_t)
    
    # Return cDSD
    if is_normalized == False:
        return X_i
    
    # Normalize with steady state
    SS = np.sqrt(d * e)
    SS = np.where(SS != 0, 1/SS, 0)
    #  SS = compute_pinverse_diagonal(np.diag(SS.flatten(

    #return np.matmul(X_i, SS)
    return X_i * SS.T

def compute_dist(X, metric = "euclidean"):
    """
    Given the embedding, compute the pairwise distance between rows.
    Default metric = "L2" or "euclidean"
    """
    return squareform(pdist(X, metric = metric))
