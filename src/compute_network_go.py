#!/cluster/tufts/cowenlab/.envs/denoise/bin/python
import sys
import os

import argparse
import numpy as np
import json
import pandas as pd
import time
import sys
import mygene
from goatools.base import get_godag
from goatools.gosubdag.gosubdag import GoSubDag
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.base import download_ncbi_associations
import itertools


####################### GOATOOLS START ###########################################
class GoTool:
    """
    Class used to process GO terms
    """
    def __init__(self, obo_location):
        """
        Location of obo file needed `go-basic.obo` 
        """
        self.godag = get_godag(obo_location, optional_attrs='relationship')
        
    def get_labels(self, filters = None):
        """
        Returns all the GO labels satisfying the certain filter. Filter is a dictionary with terms
        {'max_level' :  int #maximum level of GO terms allowed
        'min_level'  :  int #minimum level of GO terms allowed
        'namespace'  :  Namespace of the GO terms, P, F or C.
        """
        if filters == None:
            return list(self.godag.keys())
        go_terms   = []
        for k in self.godag.keys():
            k_val = self.godag[k]
            if "max_level" in filters:
                if k_val.level > filters["max_level"]:
                    continue
            if "min_level" in filters:
                if k_val.level < filters["min_level"]:
                    continue
            if "namespace" in filters:
                if k_val.namespace != filters["namespace"]:
                    continue
            go_terms.append(k)
        return go_terms

    
def get_go_labels(filter_protein, 
                  filter_label, 
                  entrez_labels, 
                  gene_to_go_file,
                  obo_file,
                  species_id,
                  anno_map = lambda x : x,
                  verbose = True):
    """
    Given a list of proteins indexed by their entrez ids, returns the GO terms involved, 
    that satisfies certain conditions outlined by filter_protein and filter_label
    
    """
    
    def log(strng):
        if verbose:
            print(strng)
            
    # Read the gene to go file, for a certain species id
    objanno = Gene2GoReader(gene_to_go_file, taxids=[species_id]) # 9606 for human
    
    # go2geneids essentially is a map that maps GO -> [genes]
    go2geneids = objanno.get_id2gos(namespace=filter_protein["namespace"], 
                                    go2geneids=True)
    
    mg = mygene.MyGeneInfo()
    
    # Use the GoTool described above, to get the filtered labels
    gt          = GoTool(obo_file)
    
    # This is the complete list of GO Labels filtered out only based on thier depth and namespace (P, F or C)
    labels      = gt.get_labels(filter_label)
    labels_dict = {}
    f_labels    = []
    
    for key in labels:
        """
        Check only the filtered labels
        """
        if key not in go2geneids:
            continue
         
        assoc_genes   = go2geneids[key]
        
        # `f_assoc_genes` to filter genes according to the filter_protein
        f_assoc_genes = list(set(assoc_genes).intersection(set(entrez_labels))) 
        
        # Removes the GO terms if it very sparsely annotates the protein list 
        lower_bound_satisfied = ("lower_bound" in filter_protein and len(f_assoc_genes) > filter_protein["lower_bound"])
        target_go_satisfied   = ("target_gos" in filter_protein and key in filter_protein["target_gos"])
        
        if lower_bound_satisfied or target_go_satisfied:
            labels_dict[key] = f_assoc_genes
            f_labels.append(key)
    log(f"Number of GO-terms: {len(f_labels)}")
    return f_labels, labels_dict
################### GOATOOLS CODE #####################################

def convert_to_entrez(genes, entrezfile = "entrezfile.tmp.json"):
    """
    Convert to entrez
    """
    if os.path.exists(entrezfile):
        with open(entrezfile, "r") as ef:
            emap = json.load(ef)
            entrez = [int(ent) for ent in emap.keys()]
            return entrez, emap

    
    mg     = mygene.MyGeneInfo()

    # compute batch of fifty
    no_genes = len(genes)
    BATCH_SIZE = 100
    no_batch = int(no_genes / BATCH_SIZE)
    print(f"Number of batches: {no_batch}...")
    gene_batch = [genes[i * BATCH_SIZE: (i + 1) * BATCH_SIZE] for i in range(no_batch)]
    gene_batch.append(genes[no_batch * BATCH_SIZE:])

    b_output = []
    for i, batch in enumerate(gene_batch):
        print(f"Processing batch {i * BATCH_SIZE}...{(i+1)*BATCH_SIZE}.")
        b_op = mg.querymany(batch, scopes = "symbol", species = 9606)
        b_output.append(b_op)
    
    output   = list(itertools.chain(*b_output))
    
    output = [op for op in output if op is not None]
    entrezlist = [int(op["entrezgene"]) for op in output
                  if "entrezgene" in op]
    print("Gene-symbol to entrez conversion complete! Saving...")
    emap = {
            op["entrezgene"]: op["_id"]
            for op in output if "entrezgene" in op
        }
    
    with open(entrezfile, "w") as of:
        json.dump(emap, of)
    
    return entrezlist, emap


def get_go_lab(go_type,
               min_level,
               min_prot,
               org_id = 9606,
               go_folder = "",
               entrez_proteins = []):
    """
    Function to get GO labels
    """
    GOT = "biological_process" if go_type == "P" else "molecular_function"
    GOT = "cellular_component" if go_type == "C" else GOT

    filter_label = {"namespace": GOT, "min_level": min_level}
    filter_prot  = {"namespace": GOT, "lower_bound": min_prot}
    labels, go_prots = get_go_labels(filter_prot,
                                     filter_label,
                                     entrez_proteins,
                                     f"../data/go/go/gene2go",
                                     f"../data/go/go/go-basic.obo",
                                     org_id,
                                     verbose = True)
    return labels, go_prots

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--network", required = True)
    parser.add_argument("--output", required = True)
    return parser.parse_args()

def main(args):
    """
    Main code
    """
    df = pd.read_csv(args.network,
                     delim_whitespace = True,
                     header = None)

    genelist     = list(set(df[0]).union(set(df[1])))
    eprots, emap = convert_to_entrez(genelist)

    mf_labels, go_mf = get_go_lab("F",
                                  5,
                                  10,
                                  entrez_proteins = eprots)

    bp_labels, go_bp = get_go_lab("P",
                                  5,
                                  10,
                                  entrez_proteins = eprots)

    cc_labels, go_cc = get_go_lab("C",
                                  5,
                                  10,
                                  entrez_proteins = eprots)

    gomaps = []
    for label in go_mf:
        for prot in [str(pt) for pt in go_mf[label] if str(pt) in emap]:
            
            gomaps.append((label, emap[prot], "molecular_function"))
            
    for label in go_bp:
        for prot in [str(pt) for pt in go_bp[label] if str(pt) in emap]:
            gomaps.append((label, emap[prot], "biological_process"))
        
    for label in go_cc:
        for prot in [str(pt) for pt in go_cc[label] if str(pt) in emap]:
            gomaps.append((label, emap[prot], "cellular_component"))

    
    df_go = pd.DataFrame(gomaps, columns = ["GO-LABEL",
                                           "PROTEIN",
                                           "GO-LABEL-TYPE"])
    df_go.to_csv(args.output, sep = "\t", index = None)
    return
            
if __name__ == "__main__":
    main(get_args())
