from dsd_utils import *
import argparse
import pandas as pd
import numpy as np

def get_args():
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


if __name__ == "__main__":
    main(get_args())