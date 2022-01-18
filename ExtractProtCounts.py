#!/usr/bin/env python3

import os
import sys

def ExtractProtCounts(filename):
    """Extract the total number of proteins within a given phylogenetic group's proteomes"""
    os.chdir("./SP_types_counts_phylogenetic_groups")

    #Open file
    try:
        infile = open(filename, "r")
    except IOError as err:
        print("Error: ", str(err))
        sys.exit(1)

    #Protein count will always be on 2nd line on which it is the 2nd tab-separated element
    infile.readline()
    data_line = infile.readline()
    prot_count = data_line.split("\t")[1]
    proteome_count = data_line.split("\t")[2]

    infile.close()
    os.chdir("../")

    return prot_count, proteome_count