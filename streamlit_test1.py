#!/usr/bin/env python3

import streamlit as st
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import math
import sys


###FUNCTIONS###



def MakeHistograms(filename, sp_type, group):
    #Test it out for single organisms
    os.chdir("./SP_regions_counts_phylogenetic_groups")

    try:
        infile = open(filename, "r")
    except IOError as err:
        print("Error: ", str(err))
        sys.exit(1)

    #Initialize
    n_region_sp = []
    h_region_sp = []
    c_region_sp = []
    n_region_lipo = []
    h_region_lipo = []
    n_region_tat = []
    h_region_tat = []
    c_region_tat = []
    n_region_tatlipo = []
    h_region_tatlipo = []
    n_region_pilin = []

    count_wrong = 0
    count_right = 0

    for line in infile:
        #Collect Sec SPI region lengths
        if line.split()[-1] == "SP":
            if line.split()[3] != "0" and line.split()[2] != "0" and line.split()[1] != "0":
                n_region_sp.append(int(line.split()[1]))
                h_region_sp.append(int(line.split()[2]))
                c_region_sp.append(int(line.split()[3]))
                count_right += 1
            else:
                #print("The signal peptide of ", line.split()[0], " has been wrongly predicted and will not be used in the analyses. Please run this analysis again.")
                count_wrong += 1
        
        #Collect Sec SPII region lengths
        elif line.split()[-1] == "LIPO":
            #Marker for throwing away wrongly classified SPs
            if line.split()[3] == "0" and line.split()[2] != "0" and line.split()[1] != "0":
                n_region_lipo.append(int(line.split()[1]))
                h_region_lipo.append(int(line.split()[2]))
                count_right += 1
            else:
                #print("The signal peptide of ", line.split()[0], " has been wrongly predicted and will not be used in the analyses. Please run this analysis again.")
                count_wrong += 1
        
        #Collect Tat SPI region lengths
        elif line.split()[-1] == "TAT":
            if line.split()[3] != "0" and line.split()[2] != "0" and line.split()[1] != "0":
                n_region_tat.append(int(line.split()[1]))
                h_region_tat.append(int(line.split()[2]))
                c_region_tat.append(int(line.split()[3]))
                count_right += 1
            else:
                #print("The signal peptide of ", line.split()[0], " has been wrongly predicted and will not be used in the analyses. Please run this analysis again.")
                count_wrong += 1
        
        #Collect Tat SPII region lengths
        elif line.split()[-1] == "TATLIPO":
            #Marker for throwing away wrongly classified SPs
            if line.split()[3] == "0" and line.split()[2] != "0" and line.split()[1] != "0":
                n_region_tatlipo.append(int(line.split()[1]))
                h_region_tatlipo.append(int(line.split()[2]))
                count_right += 1
            else:
                #print("The signal peptide of ", line.split()[0], " has been wrongly predicted and will not be used in the analyses. Please run this analysis again.")
                count_wrong += 1

        #Collect Sec SPIII region lengths
        elif line.split()[-1] == "PILIN":
            #Marker for throwing away wrongly classified SPs
            if line.split()[3] == "0" and line.split()[2] == "0" and line.split[1] != "0":
                n_region_pilin.append(int(line.split()[1]))
                count_right += 1
            else:
                #print("The signal peptide of ", line.split()[0], " has been wrongly predicted and will not be used in the analyses. Please run this analysis again.")
                count_wrong += 1

    #This can eventually be evaluated; just a checker
    #print(count_right)
    #print(count_wrong)


    def __Histograms(region, sp_type, group, region_lengths, w, hist_col):
        region_lengths.sort()
        if len(region_lengths) != 0: 
            n = math.ceil((region_lengths[-1] - region_lengths[0])/w)
            plt.style.use('ggplot')
            fig, ax = plt.subplots()
            ax.hist(region_lengths, bins = n, density = True, color = hist_col)
            ax.set_title("Distribution of "+ sp_type + " " + region + "-region lengths \n for taxonomical group " + group)
            #plt.show()
            st.pyplot(fig)

    if sp_type == "SP":
        if len(n_region_sp) == 0:
            st.write("No signal peptides of the chosen type where found for this group.")
        else:
            __Histograms("n", "sec SPI", group, n_region_sp, 1, "blue")
            __Histograms("h", "sec SPI", group, h_region_sp, 1, "green")
            __Histograms("c", "sec SPI", group, c_region_sp, 1, "red")
    elif sp_type == "LIPO":
        if len(n_region_lipo) == 0:
            st.write("No signal peptides of the chosen type where found for this group.")
        else:
            __Histograms("n", "sec SPII", group, n_region_lipo, 1, "blue")
            __Histograms("h", "sec SPII", group, h_region_lipo, 1, "green")
    elif sp_type == "TAT":
        if len(n_region_tat) == 0:
            st.write("No signal peptides of the chosen type where found for this group.")
        else:
            __Histograms("n", "tat SPI", group, n_region_tat, 1, "blue")
            __Histograms("h", "tat SPI", group, h_region_tat, 1, "green")
            __Histograms("c", "tat SPI", group, c_region_tat, 1, "red")
    elif sp_type == "TATLIPO":
        if len(n_region_tatlipo) == 0:
            st.write("No signal peptides of the chosen type where found for this group.")
        else:
            __Histograms("n", "tat SPII", group, n_region_tatlipo, 1, "blue")
            __Histograms("h", "tat SPII", group, h_region_tatlipo, 1, "green")
    elif sp_type == "PILIN":
        if len(n_region_pilin) == 0:
            st.write("No signal peptides of the chosen type where found for this group.")
        else:
            __Histograms("n", "sec SPIII", group, n_region_pilin, 1, "blue")

    infile.close()
    os.chdir("../")


def ExtractTaxonomyData():
    os.chdir("./TaxonomyData")
    files = os.listdir()

    tax_data = dict()
    all_taxa = list()
    all_taxa.append("")

    for i in range(len(files)):
        infile = open(files[i], "r")
        for line in infile:
            taxa = line.split("\t")
            tax_data[taxa[0]] = [taxa[-1], files[i][:-3]]

            all_taxa.append(taxa[0])

            if files[i] != "no rank.txt":
                all_taxa.append(taxa[-1])

        infile.close()

    os.chdir("../")

    all_taxa = tuple(all_taxa)

    return all_taxa, tax_data


all_taxa = ExtractTaxonomyData()[0]
tax_data = ExtractTaxonomyData()[1]


#########################################################################################################################################
#Streamlit design

st.title("Signal peptides across the Tree of Life")

#execute = False

option = st.selectbox(
     'Select the phylogenetic group you wish to view analyses for',
     all_taxa)
#if st.button('Click here, if ' + option + 'is the group you wish to see analyses for'):
#    execute = True
with st.expander("Further information"):
     st.write("""
         In the search box above you can pick any given taxonomical group on any given rank, 
         by searching for either the scientific name of that group or by searching for the 
         taxonomical ID as defined by NCBI's Taxonomy Database, which can be accessed here:
         https://www.ncbi.nlm.nih.gov/taxonomy
     """)

for key, value in tax_data.items():
    if option == value[0]:
        scient_name = option
        tax_id = key
    elif option == key:
        tax_id = option
        scient_name = value[0]

#full_desc = scient_name + " (" + tax_id + ")"

if option != "":
    #if execute == True:
    st.header('You are viewing the analyses for the phylogenetic group ' + scient_name)

    select_sp = st.radio("Please select the signal peptide type you wish to explore", ('Sec SPI', "Sec SPII", "Tat SPI", "Tat SPII", "Sec SPIII"))
    with st.expander("Further information"):
     st.write("""
         Write briefly about the SP types
     """)
    if select_sp == "Sec SPI":
        MakeHistograms(tax_id+"_SP_regions_counts.tab", "SP", group = scient_name)
    elif select_sp == "Sec SPII":
        MakeHistograms(tax_id+"_SP_regions_counts.tab", "LIPO", group = scient_name)
    elif select_sp == "Sec SPIII":
        MakeHistograms(tax_id+"_SP_regions_counts.tab", "PILIN", group = scient_name)
    elif select_sp == "Tat SPI":
        MakeHistograms(tax_id+"_SP_regions_counts.tab", "TAT", group = scient_name)
    elif select_sp == "Tat SPII":
        MakeHistograms(tax_id+"_SP_regions_counts.tab", "TATLIPO", group = scient_name)


