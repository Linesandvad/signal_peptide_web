#!/usr/bin/env python3

import streamlit as st
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import math
import sys


###FUNCTIONS###

def CountSPs(filename):
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

    infile.close()
    os.chdir("../")

    return n_region_sp, h_region_sp, c_region_sp, n_region_lipo, h_region_lipo, n_region_tat, h_region_tat, c_region_tat, n_region_tatlipo, h_region_tatlipo, n_region_pilin


def MakeHistograms(filename, sp_type, group):
    n_region_sp = CountSPs(filename)[0]
    h_region_sp = CountSPs(filename)[1]
    c_region_sp = CountSPs(filename)[2]
    n_region_lipo = CountSPs(filename)[3]
    h_region_lipo = CountSPs(filename)[4]
    n_region_tat = CountSPs(filename)[5]
    h_region_tat = CountSPs(filename)[6]
    c_region_tat = CountSPs(filename)[7]
    n_region_tatlipo = CountSPs(filename)[8]
    h_region_tatlipo = CountSPs(filename)[9]
    n_region_pilin = CountSPs(filename)[10]
    
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
            st.write("No signal peptides of the chosen type were found for this group.")
        elif len(n_region_sp) > 0 and len(n_region_sp) <= 5:
            st.write("Only 5 or less signal peptides of the chosen type were found for this group.")
        else:
            __Histograms("n", "sec SPI", group, n_region_sp, 1, "blue")
            __Histograms("h", "sec SPI", group, h_region_sp, 1, "green")
            __Histograms("c", "sec SPI", group, c_region_sp, 1, "red")
    elif sp_type == "LIPO":
        if len(n_region_lipo) == 0:
            st.write("No signal peptides of the chosen type were found for this group.")
        elif len(n_region_lipo) > 0 and len(n_region_lipo) <= 5:
            st.write("Only 5 or less signal peptides of the chosen type were found for this group.")
        else:
            __Histograms("n", "sec SPII", group, n_region_lipo, 1, "blue")
            __Histograms("h", "sec SPII", group, h_region_lipo, 1, "green")
    elif sp_type == "TAT":
        if len(n_region_tat) == 0:
            st.write("No signal peptides of the chosen type were found for this group.")
        elif len(n_region_tat) > 0 and len(n_region_tat) <= 5:
            st.write("Only 5 or less signal peptides of the chosen type were found for this group.")
        else:
            __Histograms("n", "tat SPI", group, n_region_tat, 1, "blue")
            __Histograms("h", "tat SPI", group, h_region_tat, 1, "green")
            __Histograms("c", "tat SPI", group, c_region_tat, 1, "red")
    elif sp_type == "TATLIPO":
        if len(n_region_tatlipo) == 0:
            st.write("No signal peptides of the chosen type were found for this group.")
        elif len(n_region_tatlipo) > 0 and len(n_region_tatlipo) <= 5:
            st.write("Only 5 or less signal peptides of the chosen type were found for this group.")
        else:
            __Histograms("n", "tat SPII", group, n_region_tatlipo, 1, "blue")
            __Histograms("h", "tat SPII", group, h_region_tatlipo, 1, "green")
    elif sp_type == "PILIN":
        if len(n_region_pilin) == 0:
            st.write("No signal peptides of the chosen type were found for this group.")
        elif len(n_region_pilin) > 0 and len(n_region_pilin) <= 5:
            st.write("Only 5 or less signal peptides of the chosen type were found for this group.")
        else:
            __Histograms("n", "sec SPIII", group, n_region_pilin, 1, "blue")


def ExtractProteinCounts(filename):
    os.chdir("./SP_types_counts_phylogenetic_groups")

    try:
        infile = open(filename, "r")
    except IOError as err:
        print("Error: ", str(err))
        sys.exit(1)

    infile.readline()
    data_line = infile.readline()
    prot_count = data_line.split("\t")[1]

    infile.close()
    os.chdir("../")

    return prot_count

#print(ExtractProteinCounts("2157_SP_types_counts.tab"))

def VisualizeSPCounts(tax_id):
    n_region_sp = CountSPs(tax_id+"_SP_regions_counts.tab")[0]
    n_region_lipo = CountSPs(tax_id+"_SP_regions_counts.tab")[3]
    n_region_tat = CountSPs(tax_id+"_SP_regions_counts.tab")[5]
    n_region_tatlipo = CountSPs(tax_id+"_SP_regions_counts.tab")[8]
    n_region_pilin = CountSPs(tax_id+"_SP_regions_counts.tab")[10]


    prot_count = ExtractProteinCounts(tax_id+"_SP_types_counts.tab")

    sp, lipo, tat, tatlipo, pilin = st.columns(5)
    sp.metric("Sec SPI", len(n_region_sp))
    lipo.metric("Sec SPII", len(n_region_lipo))
    tat.metric("Tat SPI", len(n_region_tat))
    tatlipo.metric("Tat SPII", len(n_region_tatlipo))
    pilin.metric("Sec SPIII", len(n_region_pilin))
    
def VisualizeSPFrequencies(tax_id):
    n_region_sp = CountSPs(tax_id+"_SP_regions_counts.tab")[0]
    n_region_lipo = CountSPs(tax_id+"_SP_regions_counts.tab")[3]
    n_region_tat = CountSPs(tax_id+"_SP_regions_counts.tab")[5]
    n_region_tatlipo = CountSPs(tax_id+"_SP_regions_counts.tab")[8]
    n_region_pilin = CountSPs(tax_id+"_SP_regions_counts.tab")[10]

    prot_count = ExtractProteinCounts(tax_id+"_SP_types_counts.tab")

    sp, lipo, tat, tatlipo, pilin = st.columns(5)
    sp.metric("Sec SPI", str(round((len(n_region_sp)/int(prot_count)*100), 2))+" %")
    lipo.metric("Sec SPII", str(round((len(n_region_lipo)/int(prot_count)*100), 2))+" %")
    tat.metric("Tat SPI", str(round((len(n_region_tat)/int(prot_count)*100), 2))+" %")
    tatlipo.metric("Tat SPII", str(round((len(n_region_tatlipo)/int(prot_count)*100), 2))+" %")
    pilin.metric("Sec SPIII", str(round((len(n_region_pilin)/int(prot_count)*100), 2))+" %")
    with st.expander("Further information"):
     st.write("""
        A total of""", prot_count, """ proteins were found within this phylogenetic group. 
        The fractions are calculated as the count of signal peptides of a given type divided by 
        the total number of proteins.
        """)


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





#############################################################Streamlit design############################################################################


st.title("Signal peptides across the Tree of Life")


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

if option != "":
    st.header('You are viewing the analyses for the phylogenetic group ' + scient_name)
    #st.subheader('(' + tax_id + ')')

    st.subheader("Counts of each signal peptide type")
    VisualizeSPCounts(tax_id)

    st.subheader("Frequency of proteins tagged with each type of signal peptide")
    VisualizeSPFrequencies(tax_id)

    select_sp = st.radio("Please select the signal peptide type you wish to explore", ('Sec SPI', "Sec SPII", "Tat SPI", "Tat SPII", "Sec SPIII"))
    with st.expander("Further information"):
     st.write("""
         Write briefly about the SP types?
     """)
    st.subheader("Histograms of region length distributions")
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

