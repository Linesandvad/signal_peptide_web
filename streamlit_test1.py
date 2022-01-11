#!/usr/bin/env python3

import streamlit as st
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import math
import sys
import pickle

#Remember to by en main directory with all subdirectories

#######################FUNCTIONS###############################


def CountSPs(filename):
    """Count region lengths for the different SP types within a phylogenetic group"""
    os.chdir("./SP_regions_counts_phylogenetic_groups")

    #Open file with SP region length counts
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

    #Control check; remove wrongly predicted signal peptides
    count_wrong = 0
    count_right = 0

    #Go through all SPs, one at a time
    for line in infile:
        #Collect Sec SPI region lengths
        if line.split()[-1] == "SP":
            #Filter wrongly predicted SPs from analyses
            if line.split()[3] != "0" and line.split()[2] != "0" and line.split()[1] != "0":
                #Append region length observations to separate lists
                n_region_sp.append(int(line.split()[1]))
                h_region_sp.append(int(line.split()[2]))
                c_region_sp.append(int(line.split()[3]))
                count_right += 1
            else:
                #print("The signal peptide of ", line.split()[0], " has been wrongly predicted and will not be used in the analyses. Please run this analysis again.")
                count_wrong += 1
        
        #Collect Sec SPII region lengths
        elif line.split()[-1] == "LIPO":
            #Filter wrongly predicted SPs from analyses
            if line.split()[3] == "0" and line.split()[2] != "0" and line.split()[1] != "0":
                #Append region length observations to separate lists
                n_region_lipo.append(int(line.split()[1]))
                h_region_lipo.append(int(line.split()[2]))
                count_right += 1
            else:
                #print("The signal peptide of ", line.split()[0], " has been wrongly predicted and will not be used in the analyses. Please run this analysis again.")
                count_wrong += 1
        
        #Collect Tat SPI region lengths
        elif line.split()[-1] == "TAT":
            #Filter wrongly predicted SPs from analyses
            if line.split()[3] != "0" and line.split()[2] != "0" and line.split()[1] != "0":
                #Append region length observations to separate lists
                n_region_tat.append(int(line.split()[1]))
                h_region_tat.append(int(line.split()[2]))
                c_region_tat.append(int(line.split()[3]))
                count_right += 1
            else:
                #print("The signal peptide of ", line.split()[0], " has been wrongly predicted and will not be used in the analyses. Please run this analysis again.")
                count_wrong += 1
        
        #Collect Tat SPII region lengths
        elif line.split()[-1] == "TATLIPO":
            #Filter wrongly predicted SPs from analyses
            if line.split()[3] == "0" and line.split()[2] != "0" and line.split()[1] != "0":
                #Append region length observations to separate lists
                n_region_tatlipo.append(int(line.split()[1]))
                h_region_tatlipo.append(int(line.split()[2]))
                count_right += 1
            else:
                #print("The signal peptide of ", line.split()[0], " has been wrongly predicted and will not be used in the analyses. Please run this analysis again.")
                count_wrong += 1

        #Collect Sec SPIII region lengths
        elif line.split()[-1] == "PILIN":
            #Filter wrongly predicted SPs from analyses
            if line.split()[3] == "0" and line.split()[2] == "0":
                #Append region length observations to separate lists
                n_region_pilin.append(int(line.split()[1]))
                count_right += 1
            else:
                #print("The signal peptide of ", line.split()[0], " has been wrongly predicted and will not be used in the analyses. Please run this analysis again.")
                count_wrong += 1

    infile.close()
    os.chdir("../")

    #Return all region counts
    return n_region_sp, h_region_sp, c_region_sp, n_region_lipo, h_region_lipo, n_region_tat, h_region_tat, c_region_tat, n_region_tatlipo, h_region_tatlipo, n_region_pilin


def MakeHistograms(filename, sp_type, group):
    """Make histograms of SP region length distributions, given a phylogenetic group and SP type"""

    #Collect data from function CountSPs
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
        """Make histograms, display to streamlit"""
        #Sort region lengths from shortest to longest
        region_lengths.sort()
        #Only produce non-empty histograms
        if len(region_lengths) != 0: 
            #Define number of bins (1 bin for each integer)
            n = math.ceil((region_lengths[-1] - region_lengths[0])/w)
            plt.style.use('ggplot')
            fig, ax = plt.subplots()
            #Create histogram
            ax.hist(region_lengths, bins = n, density = True, color = hist_col)
            #Create histogram title
            ax.set_title("Distribution of "+ sp_type + " " + region + "-region lengths \n for taxonomical group " + group)
            #Display histogram in streamlit app
            st.pyplot(fig)

    #Use this function if all region counts in a region for a group is the same (rare case)
    def __HistogramsOneLength(region, sp_type, group, region_lengths, hist_col):
        """Make histograms, display to streamlit"""
        #Sort region lengths from shortest to longest
        region_lengths.sort()
        #Only produce non-empty histograms
        if len(region_lengths) != 0: 
            #Define number of bins (1 bin for each integer)
            n = 1
            plt.style.use('ggplot')
            fig, ax = plt.subplots()
            #Create histogram
            ax.hist(region_lengths, bins = n, density = True, color = hist_col)
            #Create histogram title
            ax.set_title("Distribution of "+ sp_type + " " + region + "-region lengths \n for taxonomical group " + group)
            #Display histogram in streamlit app
            st.pyplot(fig)

    #Display histograms of SP region lengths of the selected SP type
    #Displaying histograms for Sec SPI
    if sp_type == "SP":
        if len(n_region_sp) == 0:
            st.info("No signal peptides of the chosen type were found for this group.")
        else:
            try:
                __Histograms("n", "sec SPI", group, n_region_sp, 1, "blue")
            except ValueError as err:
                __HistogramsOneLength("n", "sec SPI", group, n_region_sp, "blue")
                if len(n_region_sp) > 1:
                    st.info("All Sec SPI n-regions within this group consist of " + str(n_region_sp[0]) + " residues.")
            try:
                __Histograms("h", "sec SPI", group, h_region_sp, 1, "green")
            except ValueError as err:
                __HistogramsOneLength("h", "sec SPI", group, h_region_sp, "green")
                if len(h_region_sp) > 1:
                    st.info("All Sec SPI h-regions within this group consist of " + str(h_region_sp[0]) + " residues.")
            try:
                __Histograms("c", "sec SPI", group, c_region_sp, 1, "red")
            except ValueError as err:
                __HistogramsOneLength("c", "sec SPI", group, c_region_sp, "red")
                if len(c_region_sp) > 1:
                    st.info("All Sec SPI c-regions within this group consist of " + str(c_region_sp[0]) + " residues.")
    
    #Displaying histograms for Sec SPII
    elif sp_type == "LIPO":
        if len(n_region_lipo) == 0:
            st.info("No signal peptides of the chosen type were found for this group.")
        else:
            try:
                __Histograms("n", "sec SPII", group, n_region_lipo, 1, "blue")
            except ValueError as err:
                __HistogramsOneLength("n", "sec SPII", group, n_region_lipo, "blue")
                if len(n_region_lipo) > 1:
                    st.info("All Sec SPII n-regions within this group consist of " + str(n_region_lipo[0]) + " residues.")
            try:
                __Histograms("h", "sec SPII", group, h_region_lipo, 1, "green")
            except ValueError as err:
                __HistogramsOneLength("h", "sec SPII", group, h_region_lipo, "green")
                if len(h_region_lipo) > 1:
                    st.info("All Sec SPII h-regions within this group consist of " + str(h_region_lipo[0]) + " residues.")
    
    #Displaying histograms for Tat SPI
    elif sp_type == "TAT":
        if len(n_region_tat) == 0:
            st.info("No signal peptides of the chosen type were found for this group.")
        else:
            try:
                __Histograms("n", "tat SPI", group, n_region_tat, 1, "blue")
            except ValueError as err:
                __HistogramsOneLength("n", "tat SPI", group, n_region_tat, "blue")
                if len(n_region_tat) > 1:
                    st.info("All Tat SPI n-regions within this group consist of " + str(n_region_tat[0]) + " residues.")
            try:
                __Histograms("h", "tat SPI", group, h_region_tat, 1, "green")
            except ValueError as err:
                __HistogramsOneLength("h", "tat SPI", group, h_region_tat, "green")
                if len(h_region_tat) > 1:
                    st.info("All Tat SPI h-regions within this group consist of " + str(h_region_tat[0]) + " residues.")
            try:
                __Histograms("c", "tat SPI", group, c_region_tat, 1, "red")
            except ValueError as err:
                __HistogramsOneLength("c", "tat SPI", group, c_region_tat, "red")
                if len(c_region_tat) > 1:
                    st.info("All Tat SPI c-regions within this group consist of " + str(c_region_tat[0]) + " residues.")
    
    #Displaying histograms for Tat SPII
    elif sp_type == "TATLIPO":
        if len(n_region_tatlipo) == 0:
            st.info("No signal peptides of the chosen type were found for this group.")
        else:
            try:
                __Histograms("n", "tat SPII", group, n_region_tatlipo, 1, "blue")
            except ValueError as err:
                __HistogramsOneLength("n", "tat SPII", group, n_region_tatlipo, "blue")
                if len(n_region_tatlipo) > 1:
                    st.info("All Tat SPII n-regions within this group consist of " + str(n_region_tatlipo[0]) + " residues.")
            try:
                __Histograms("h", "tat SPII", group, h_region_tatlipo, 1, "green")
            except ValueError as err:
                __HistogramsOneLength("h", "tat SPII", group, h_region_tatlipo, "green")
                if len(h_region_tatlipo) > 1:
                    st.info("All Tat SPII h-regions within this group consist of " + str(h_region_tatlipo[0]) + " residues.")
    
    #Displaying histograms for Sec SPIII
    elif sp_type == "PILIN":
        if len(n_region_pilin) == 0:
            st.info("No signal peptides of the chosen type were found for this group.")
        else:
            try:
                __Histograms("n", "sec SPIII", group, n_region_pilin, 1, "blue")
            except ValueError as err:
                __HistogramsOneLength("n", "sec SPIII", group, n_region_pilin, "blue")
                if len(n_region_pilin) > 1:
                    st.info("All Sec SPIII signal peptides within this group contains " + str(n_region_pilin[0]) + " residues.")


def ExtractProteinCounts(filename):
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


def VisualizeProtCounts(tax_id):
    """Display protein- and proteome counts for a given taxonomical group"""
    #Save protein and proteome counts
    prot_count = ExtractProteinCounts(tax_id+"_SP_types_counts.tab")[0]
    proteome_count = ExtractProteinCounts(tax_id+"_SP_types_counts.tab")[1]

    proteome, protein = st.columns(2)

    proteome.metric("Proteome count", proteome_count) 
    protein.metric("Protein count", prot_count) 



def VisualizeSPCounts(tax_id):
    """Displays the counts of each SP type for a given group to streamlit"""
    #Save region counts
    n_region_sp = CountSPs(tax_id+"_SP_regions_counts.tab")[0]
    n_region_lipo = CountSPs(tax_id+"_SP_regions_counts.tab")[3]
    n_region_tat = CountSPs(tax_id+"_SP_regions_counts.tab")[5]
    n_region_tatlipo = CountSPs(tax_id+"_SP_regions_counts.tab")[8]
    n_region_pilin = CountSPs(tax_id+"_SP_regions_counts.tab")[10]

    #Setup in streamlit
    sp, lipo, tat, tatlipo, pilin = st.columns(5)
    #Display count of each SP type to streamlit
    sp.metric("Sec SPI", len(n_region_sp))
    lipo.metric("Sec SPII", len(n_region_lipo))
    tat.metric("Tat SPI", len(n_region_tat))
    tatlipo.metric("Tat SPII", len(n_region_tatlipo))
    pilin.metric("Sec SPIII", len(n_region_pilin))
    
def VisualizeSPFrequencies(tax_id):
    """Displays the fractions of proteins tagged with each SP type for a given group to streamlit"""
    #Save region counts
    n_region_sp = CountSPs(tax_id+"_SP_regions_counts.tab")[0]
    n_region_lipo = CountSPs(tax_id+"_SP_regions_counts.tab")[3]
    n_region_tat = CountSPs(tax_id+"_SP_regions_counts.tab")[5]
    n_region_tatlipo = CountSPs(tax_id+"_SP_regions_counts.tab")[8]
    n_region_pilin = CountSPs(tax_id+"_SP_regions_counts.tab")[10]

    #Save total protein count for phylogenetic group
    prot_count = ExtractProteinCounts(tax_id+"_SP_types_counts.tab")[0]
    
    try:
        #Setup in streamlit
        sp, lipo, tat, tatlipo, pilin = st.columns(5)
        #Display fractions of each SP type to streamlit, round to 2 significant digits
        sp.metric("Sec SPI", str(round((len(n_region_sp)/int(prot_count)*100), 2))+" %")
        lipo.metric("Sec SPII", str(round((len(n_region_lipo)/int(prot_count)*100), 2))+" %")
        tat.metric("Tat SPI", str(round((len(n_region_tat)/int(prot_count)*100), 2))+" %")
        tatlipo.metric("Tat SPII", str(round((len(n_region_tatlipo)/int(prot_count)*100), 2))+" %")
        pilin.metric("Sec SPIII", str(round((len(n_region_pilin)/int(prot_count)*100), 2))+" %")
    except ZeroDivisionError as err:
        st.write("0 proteins where detected.")

    #Write explanation to calculations to streamlit
    with st.expander("Further information"):
     st.write("""
        The fractions of proteins tagged with a given signal peptide type
        is calculated as the count of the respective signal peptide type, divided
        by the total number of proteins in this phylogenetic group. The total number 
        of proteins is calculated based on the number of proteins found in the proteomes 
        belonging to this phylogenetic group. The total number of
        proteins is """, prot_count, """ in this group.""")

def TreeSearchDict(domain):
    infile_dict = open(domain+"_tax_groups_children.pk", "rb")
    tree_search_file = pickle.load(infile_dict)

    return tree_search_file

tree_search_file = TreeSearchDict("archaea")

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
            #print(taxa)
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


all_taxa = ExtractTaxonomyData()[0]
tax_data = ExtractTaxonomyData()[1]



#############################################################Streamlit design############################################################################


st.title("Signal peptides across the Tree of Life")

option = st.selectbox(
     'Select the phylogenetic group you wish to view analyses for',
     all_taxa)
with st.expander("Further information"):
     st.write("""
         In the search box above you can pick any given taxonomical group on any given rank, 
         by searching for either the scientific name of that group or by searching for the 
         taxonomical ID as defined by NCBI's Taxonomy Database, which can be accessed here:
         https://www.ncbi.nlm.nih.gov/taxonomy
     """)


for key, value in tax_data.items():
    if option in value[1:]:
        name = option
        tax_id = key
        #Define scientific name for tree search
        for tax_name in value:
            if tax_name in tree_search_file.keys():
                scient_name = tax_name
    elif option == key:
        tax_id = option
        name = value[1]
        #Define scientific name for tree search
        for tax_name in value:
            if tax_name in tree_search_file.keys():
                scient_name = tax_name

counter = 0
if option != "":
    #Make Tree-like structure in sidebar if possible
    #Make dict where each Tax ID is key with its' children as values for Tree-like structure in sidebar
    try:
        children = tree_search_file[scient_name]
        children[0] = children[0] + " ("+ scient_name + ")"
        counter += 1
        st.sidebar.header("Tree-sctructured search for undergroups")
        option = st.sidebar.radio(
        "If you would like to look into one of the groups belonging to "+ scient_name + ", click on this groups' name.", children,
        key = counter)
        with st.sidebar.expander("Further information"):
            st.write("""
            Here you can dig further into one of the groups belonging to """ + scient_name + """ if you wish so. When clicking 
            on one of these groups, the sidebar will change and direct you into this undergroup, showing analyses for this group
            and letting you explore further undergroups if you wish to. 
         """)
        
        while option != children[0] and len(children) > 1:
            for key, value in tax_data.items():
                if option in value[1:]:
                    name = option
                    tax_id = key
                    #Define scientific name for tree search
                    for tax_name in value:
                        if tax_name in tree_search_file.keys():
                            scient_name = tax_name
                elif option == key:
                    tax_id = option
                    name = value[1]
                    #Define scientific name for tree search
                    for tax_name in value:
                        if tax_name in tree_search_file.keys():
                            scient_name = tax_name
            children = tree_search_file[scient_name]
            children[0] = children[0] + " ("+ scient_name + ")"
            counter += 1
            st.sidebar.subheader("Currently viewing "+ scient_name)
            option = st.sidebar.radio("Please click on an undergroup, to see the analyses of this group.", children,
            key = counter)
            for key, value in tax_data.items():
                if option in value[1:]:
                    name = option
                    tax_id = key
                    #Define scientific name for tree search
                    for tax_name in value:
                        if tax_name in tree_search_file.keys():
                            scient_name = tax_name
                elif option == key:
                    tax_id = option
                    name = value[1]
                    #Define scientific name for tree search
                    for tax_name in value:
                        if tax_name in tree_search_file.keys():
                            scient_name = tax_name

        st.header('You are viewing the analyses for the phylogenetic group ' + scient_name + ' ('+ tax_id + ')')

        st.subheader("Proteome and protein counts")
        VisualizeProtCounts(tax_id)

        st.subheader("Counts of each signal peptide type")
        VisualizeSPCounts(tax_id)

        st.subheader("Fraction of proteins tagged with each signal peptide type")
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


    except NameError as err:
        st.warning('The program was not able to evoke analyses of ' + option + '. Please try with another taxonomical group or organism name.')


