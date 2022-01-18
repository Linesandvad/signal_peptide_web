#!/usr/bin/env python3

import os
import sys
from ExtractProtCounts import *

def ExtractTaxonomyData():
    os.chdir("./TaxonomyData")
    #List of all taxonomic ranks
    files_tax_ranks = os.listdir()

    #Initialize 
    tax_data = dict()
    all_taxa = list()
    all_tax_ids = list()

    #Append empty object to initialize webpage without a selected group
    all_taxa.append("")

    #Loop over each filename (one filename pr. taxonomic rank)
    for i in range(len(files_tax_ranks)):
        #Open file with information about taxonomical group at the given rank
        infile = open(files_tax_ranks[i], "r")
        #Each line in file corresponds to the tax ID and name of a taxonomical group on that rank
        for line in infile:
            taxa = line.split("\t")             #Split line into Tax ID and name
            os.chdir("../")
            proteome_count = ExtractProtCounts(taxa[0]+"_SP_types_counts.tab")[1]
            os.chdir("./TaxonomyData")
            if int(proteome_count) > 0:
                #Add information to dict about taxonomical groups
                #This adds; tax ID as key, list of the tax group scientific name and tax rank as value
                tax_data[taxa[0]] = [files_tax_ranks[i][:-4]]
                if len(taxa) > 2:
                    for j in range(len(taxa)-2):
                        tax_data[taxa[0]].append(taxa[j+1])
                
                elif len(taxa) == 2:
                    tax_data[taxa[0]].append(taxa[-1][:-1])

                #append all tax IDs to list
                all_taxa.append(taxa[0])
                all_tax_ids.append(taxa[0])

                #append all tax group names to list, except the ones without rank
                if files_tax_ranks[i] != "no rank.txt":
                    if len(taxa) > 2:
                        for j in range(len(taxa)-2):
                            all_taxa.append(taxa[j+1])

                    elif len(taxa) == 2:
                        all_taxa.append(taxa[-1][:-1])

        infile.close()


    all_taxa = tuple(all_taxa)
    os.chdir("../")
    return all_taxa, tax_data