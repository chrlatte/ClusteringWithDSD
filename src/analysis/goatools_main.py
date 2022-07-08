import numpy as np
import pandas as pd
# from matplotlib.pyplot import rc_context



from goatools.base import download_go_basic_obo
from goatools.base import download_ncbi_associations
from goatools.obo_parser import GODag
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS

from genes_ncbi_test_proteincoding import GENEID2NT as GeneID2nt_mus # the filename is what we made before


def main():


    obo_fname = download_go_basic_obo() # TODO what is this 
    # print(obo_fname)
    fin_gene2go = download_ncbi_associations()
    # print(fin_gene2go)
    obodag = GODag("go-basic.obo") # dictionary of different GO terms
    # print(obodag)

    # print(type(GeneID2nt_mus)) # this is a dictionary that maps a gene id to info about it

    mapper = {} # dictionary that maps a gene symbol (ie LOC124906745) to it's gene id

    for key in GeneID2nt_mus:
        mapper[GeneID2nt_mus[key].Symbol] = GeneID2nt_mus[key].GeneID

    print(mapper)

if __name__ == "__main__":
    main()