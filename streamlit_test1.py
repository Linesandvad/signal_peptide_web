#!/usr/bin/env python3

import streamlit as st
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import math
import sys
import pickle
import logomaker
from PIL import Image

##Functions required to run app
from CountSPs import *
from MakeHistograms import *
from ExtractProtCounts import *
from VisualizeProtCounts import *
from VisualizeSPCounts import *
from VisualizeSPFrequencies import *
from TreeSearchDict import *
from ExtractTaxonomyData import *
from MakeLogos import *


#Remember to be in main directory with all subdirectories

#######################Collect data from functions###############################

#Get dict of taxonomy data for tree-structured search
tree_search_file = TreeSearchDict("archaea")

#Extract list of all taxa for search field and dict og taxonomy data connections
all_taxa = ExtractTaxonomyData()[0]
tax_data = ExtractTaxonomyData()[1]


#Insert images
seq_logo_color_schemes = Image.open('seq_logo_color_schemes.PNG')

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
        clean_rank = value[0]
        #Define scientific name for tree search
        for tax_name in value:
            if tax_name in tree_search_file.keys():
                scient_name = tax_name
    elif option == key:
        tax_id = option
        name = value[1]
        clean_rank = value[0]
        #Define scientific name for tree search
        for tax_name in value:
            if tax_name in tree_search_file.keys():
                scient_name = tax_name

counter = 0
if option != "":
    #Make Tree-like structure in sidebar if possible
    #Make dict where each Tax ID is key with its' children as values for Tree-like structure in sidebar
    try:
        children_names = tree_search_file[scient_name][::2]
        children_ranks = tree_search_file[scient_name][1::2]

        children = []
        for i in range(len(children_names)):
            child_name = children_names[i] + " "+children_ranks[i]
            children.append(child_name)

        children[0] = children[0] + "("+ scient_name + ")"

        if len(children) == 2 and "(strain)" in children[1]:
            stop = True
        else:
            counter += 1
            st.sidebar.header("Tree-structured search for undergroups")
            option = st.sidebar.radio(
            "If you would like to look into one of the groups belonging to "+ scient_name + ", click on this groups' name.", children,
            key = counter)
            with st.sidebar.expander("Further information"):
                st.write("""
                Here you can dig further into one of the groups belonging to """ + scient_name + """ if you wish so. When clicking 
                on one of these groups, the sidebar will change and direct you into this undergroup, showing analyses for this group
                and letting you explore further undergroups if you wish to. 
             """)

            stop = False

        
        while option != children[0] and not stop:
            #clean_rank_proces = option.split(" (")[1]
            #clean_rank = clean_rank_proces.split(")")[0]
            option = option.split(" (")[0]

            for key, value in tax_data.items():
                if option in value[1:]:
                    name = option
                    tax_id = key
                    clean_rank = value[0]
                    #Define scientific name for tree search
                    for tax_name in value:
                        if tax_name in tree_search_file.keys():
                            scient_name = tax_name
                elif option == key:
                    tax_id = option
                    name = value[1]
                    clean_rank = value[0]
                    #Define scientific name for tree search
                    for tax_name in value:
                        if tax_name in tree_search_file.keys():
                            scient_name = tax_name
            
            children_names = tree_search_file[scient_name][::2]
            children_ranks = tree_search_file[scient_name][1::2]

            children = []
            for i in range(len(children_names)):
                child_name = children_names[i] + " "+children_ranks[i]
                children.append(child_name)
            if len(children) == 1: 
                stop = True
            if len(children) == 2 and "(strain)" in children[1]:
                stop = True
            else:
                children[0] = children[0] + " ("+ scient_name + ")"
                counter += 1
                st.sidebar.subheader("Currently viewing "+ scient_name)
                option = st.sidebar.radio("Please click on an undergroup, to see the analyses of this group.", children,
                key = counter)

                #Split string to extract name, rank and ID of selected group
                option_redo = option.split(" (")[0]

                stop = True
                for key, value in tax_data.items():
                    if option_redo in value[1:]:
                        stop = False
                        name = option_redo                          #Save name (selected) of group
                        tax_id = key                                #Save tax ID of the group
                        clean_rank = value[0]                       #Save rank of the group
                        #Define scientific name for tree search
                        for tax_name in value:
                            if tax_name in tree_search_file.keys():
                                scient_name = tax_name              #Save scientific name
                    elif option_redo == key:
                        stop = False
                        tax_id = option_redo                        #Save tax ID of the group
                        name = value[1]                             #Save name (selected) of group
                        clean_rank = value[0]                       #Save rank of the group
                        #Define scientific name for tree search
                        for tax_name in value:
                            if tax_name in tree_search_file.keys():
                                scient_name = tax_name              #Save scientific name


        st.header('You are viewing the analyses for the ' + clean_rank + ' '+ scient_name + ' ('+ tax_id + ')')

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

        st.subheader("You have chosen to view analyses of " + select_sp)
        
        #Sec SPI analyses
        if select_sp == "Sec SPI":
            select_analyses = st.radio("Which analyses would you like to take a look at?", ("Full signal peptide length distribution", "Length distributions of each signal peptide region", "Sequence logos"))
            if select_analyses == "Full signal peptide length distribution":
                st.subheader("Histogram of signal peptide length distribution")
                MakeHistograms(tax_id+"_SP_regions_counts.tab", "SP", group = scient_name, see_regions = False, rank = clean_rank)
            elif select_analyses == "Length distributions of each signal peptide region":
                st.subheader("Histograms of region length distributions")
                MakeHistograms(tax_id+"_SP_regions_counts.tab", "SP", group = scient_name, see_regions = True, rank = clean_rank)
            elif select_analyses == "Sequence logos":
                st.subheader("Sequence logos aligned at region borders")
                
                option1,option2 = st.columns(2)
                option1_inf,option2_inf = st.columns(2)
                color_scheme = option1.radio("Select a color scheme for the sequence logos", ("Default", "Hydrophobicity", "Charge", "Chemistry"))
                with option1_inf.expander("Information about color schemes"):
                    st.image(seq_logo_color_schemes)
                    st.write("""
                        Default color scheme divides amino acids into negatively charged (red), 
                        positively charged (blue), hydrophobic (black) and others (green). This color scheme
                        is useful when examining important signal peptides properties such as
                        the positive charge of n-regions and hydrophobicity of h-regions. 
                        """)

                keep_zero = option2.radio("Do you want to include 0 as the position before alignment site?", ("Yes", "No"))
                with option2_inf.expander("Information about position assignments"):
                    st.write("""
                        Choosing either of these options do not influence the logos themselves, 
                        except changing whether the position before alignment site should be 
                        written as position 0 or position -1. 
                        """)
                #st.caption("Alignment at N-H region border")
                min_pos_n = ExtractAlignmentPosition(tax_id, "SP", "N")[0]
                len_seq_n = ExtractAlignmentPosition(tax_id, "SP", "N")[1]
                max_pos_n = int(len_seq_n) - int(min_pos_n)
                st.write("Select the range of positions in the logo aligned at the n-h region border")
                low,high = st.columns(2)
                lower_n = low.number_input('Lower value', min_value = -int(min_pos_n), max_value = -1, value = -5, step = 1)
                upper_n = high.number_input('Upper value', min_value = 1, max_value = max_pos_n, value = 7, step = 1)
                MakeLogo(tax_id,"SP","n", lower_n, upper_n, scient_name, color_scheme, keep_zero)
                #except ValueError:
                #    st.warning("Your selected value for lower part of the range is too small. Please try with a greater value.")
                min_pos_h = ExtractAlignmentPosition(tax_id, "SP", "H")[0]
                len_seq_h = ExtractAlignmentPosition(tax_id, "SP", "H")[1]
                max_pos_h = int(len_seq_h) - int(min_pos_h)
                st.write("Select the range of positions in the logo aligned at the h-c region border")
                low,high = st.columns(2)
                lower_h = low.number_input('Lower value', min_value = -int(min_pos_h), max_value = -1, value = -15, step = 1)
                upper_h = high.number_input('Upper value', min_value = 1, max_value = max_pos_h, value = 5, step = 1)
                MakeLogo(tax_id,"SP","h", lower_h, upper_h, scient_name, color_scheme, keep_zero)

                min_pos_c = ExtractAlignmentPosition(tax_id, "SP", "C")[0]
                len_seq_c = ExtractAlignmentPosition(tax_id, "SP", "C")[1]
                max_pos_c = int(len_seq_c) - int(min_pos_c)
                st.write("Select the range of positions in the logo aligned at cleavage site")
                low,high = st.columns(2)
                lower_c = low.number_input('Lower value', min_value = -int(min_pos_c), max_value = -1, value = -7, step = 1)
                upper_c = high.number_input('Upper value', min_value = 1, max_value = max_pos_c, value = 5, step = 1)
                MakeLogo(tax_id,"SP","c", lower_c, upper_c, scient_name, color_scheme, keep_zero)

                with st.expander("Creation of the sequence logos"):
                    st.write("""
                        Write something about how the logos are created (Shannon, Schneider and Steven, 1990).
                        """)
                with st.expander("Calculation of errorbars"):
                    st.write("""
                        The magnitude of the errorbars are based on the small sample correction
                        factor. The small sample correction factor for sequence logos with 
                        amino acids is calculated as
                        """)
                    st.latex(r'''e_n = \frac{1}{\ln{2}} \cdot \frac{20-1}{2\cdot n}''')
                    st.write("""
                        The height of the errorbars correspond to twice the small sample correction
                        factor. 
                        """)


                        
        #Sec SPII analyses
        elif select_sp == "Sec SPII":
            select_analyses = st.radio("Which analyses would you like to take a look at?", ("Full signal peptide length distribution", "Length distributions of each signal peptide region", "Sequence logos"))
            if select_analyses == "Full signal peptide length distribution":
                st.subheader("Histogram of signal peptide length distribution")
                MakeHistograms(tax_id+"_SP_regions_counts.tab", "LIPO", group = scient_name, see_regions = False, rank = clean_rank)
            elif select_analyses == "Length distributions of each signal peptide region":
                st.subheader("Histograms of region length distributions")
                MakeHistograms(tax_id+"_SP_regions_counts.tab", "LIPO", group = scient_name, see_regions = True, rank = clean_rank)
            elif select_analyses == "Sequence logos":
                st.subheader("Sequence logos aligned at region borders")
                option1,option2 = st.columns(2)
                option1_inf,option2_inf = st.columns(2)
                color_scheme = option1.radio("Select a color scheme for the sequence logos", ("Default", "Hydrophobicity", "Charge", "Chemistry"))
                with option1_inf.expander("Information about color schemes"):
                    st.image(seq_logo_color_schemes)
                    st.write("""
                        Default color scheme divides amino acids into negatively charged (red), 
                        positively charged (blue), hydrophobic (black) and others (green). This color scheme
                        is useful when examining important signal peptides properties such as
                        the positive charge of n-regions and hydrophobicity of h-regions. 
                        """)

                keep_zero = option2.radio("Do you want to include 0 as the position before alignment site?", ("Yes", "No"))
                with option2_inf.expander("Information about position assignments"):
                    st.write("""
                        Choosing either of these options do not influence the logos themselves, 
                        except changing whether the position before alignment site should be 
                        written as position 0 or position -1. 
                        """)
                #st.caption("Alignment at N-H region border")
                min_pos_n = ExtractAlignmentPosition(tax_id, "LIPO", "N")[0]
                len_seq_n = ExtractAlignmentPosition(tax_id, "LIPO", "N")[1]
                max_pos_n = int(len_seq_n) - int(min_pos_n)
                st.write("Select the range of positions in the logo aligned at the n-h region border")
                low,high = st.columns(2)
                lower_n = low.number_input('Lower value', min_value = -int(min_pos_n), max_value = -1, value = -5, step = 1)
                upper_n = high.number_input('Upper value', min_value = 1, max_value = max_pos_n, value = 7, step = 1)
                MakeLogo(tax_id,"LIPO","n", lower_n, upper_n, scient_name, color_scheme, keep_zero)
                #except ValueError:
                #    st.warning("Your selected value for lower part of the range is too small. Please try with a greater value.")
                min_pos_h = ExtractAlignmentPosition(tax_id, "LIPO", "H")[0]
                len_seq_h = ExtractAlignmentPosition(tax_id, "LIPO", "H")[1]
                max_pos_h = int(len_seq_h) - int(min_pos_h)
                st.write("Select the range of positions in the logo aligned at the h-c region border")
                low,high = st.columns(2)
                lower_h = low.number_input('Lower value', min_value = -int(min_pos_h), max_value = -1, value = -7, step = 1)
                upper_h = high.number_input('Upper value', min_value = 1, max_value = max_pos_h, value = 5, step = 1)
                MakeLogo(tax_id,"LIPO","h", lower_h, upper_h, scient_name, color_scheme, keep_zero)

                with st.expander("Creation of the sequence logos"):
                    st.write("""
                        Write something about how the logos are created (Shannon, Schneider and Steven, 1990).
                        """)
                with st.expander("Calculation of errorbars"):
                    st.write("""
                        The magnitude of the errorbars are based on the small sample correction
                        factor. The small sample correction factor for sequence logos with 
                        amino acids is calculated as
                        """)
                    st.latex(r'''e_n = \frac{1}{\ln{2}} \cdot \frac{20-1}{2\cdot n}''')
                    st.write("""
                        The height of the errorbars correspond to twice the small sample correction
                        factor. 
                        """)
        
        #Sec SPIII analyses
        elif select_sp == "Sec SPIII":
            select_analyses = st.radio("Which analyses would you like to take a look at?", ("Full signal peptide length distribution", "Sequence logo"))
            with st.expander("Further information"):
                st.write("""
                    The Sec SPIII signal peptide only contains one region, which is why histograms of region length distributions cannot be viewed for this signal peptide type.
                    """)

            if select_analyses == "Full signal peptide length distribution":
                st.subheader("Histogram of signal peptide length distribution")
                MakeHistograms(tax_id+"_SP_regions_counts.tab", "PILIN", group = scient_name, see_regions = False, rank = clean_rank)
            elif select_analyses == "Sequence logo":
                st.subheader("Sequence logos aligned at region borders")
                option1,option2 = st.columns(2)
                option1_inf,option2_inf = st.columns(2)
                color_scheme = option1.radio("Select a color scheme for the sequence logos", ("Default", "Hydrophobicity", "Charge", "Chemistry"))
                with option1_inf.expander("Information about color schemes"):
                    st.image(seq_logo_color_schemes)
                    st.write("""
                        Default color scheme divides amino acids into negatively charged (red), 
                        positively charged (blue), hydrophobic (black) and others (green). This color scheme
                        is useful when examining important signal peptides properties such as
                        the positive charge of n-regions and hydrophobicity of h-regions. 
                        """)

                keep_zero = option2.radio("Do you want to include 0 as the position before alignment site?", ("Yes", "No"))
                with option2_inf.expander("Information about position assignments"):
                    st.write("""
                        Choosing either of these options do not influence the logos themselves, 
                        except changing whether the position before alignment site should be 
                        written as position 0 or position -1. 
                        """)
                #st.caption("Alignment at N-H region border")
                min_pos_n = ExtractAlignmentPosition(tax_id, "PILIN", "N")[0]
                len_seq_n = ExtractAlignmentPosition(tax_id, "PILIN", "N")[1]
                max_pos_n = int(len_seq_n) - int(min_pos_n)
                st.write("Select the range of positions in the logo aligned at cleavage site")
                low,high = st.columns(2)
                lower_n = low.number_input('Lower value', min_value = -int(min_pos_n), max_value = -1, value = -7, step = 1)
                upper_n = high.number_input('Upper value', min_value = 1, max_value = max_pos_n, value = 12, step = 1)
                MakeLogo(tax_id,"PILIN","n", lower_n, upper_n, scient_name, color_scheme, keep_zero)
                
                with st.expander("Creation of the sequence logos"):
                    st.write("""
                        Write something about how the logos are created (Shannon, Schneider and Steven, 1990).
                        """)
                with st.expander("Calculation of errorbars"):
                    st.write("""
                        The magnitude of the errorbars are based on the small sample correction
                        factor. The small sample correction factor for sequence logos with 
                        amino acids is calculated as
                        """)
                    st.latex(r'''e_n = \frac{1}{\ln{2}} \cdot \frac{20-1}{2\cdot n}''')
                    st.write("""
                        The height of the errorbars correspond to twice the small sample correction
                        factor. 
                        """)

        #Tat SPI analyses
        elif select_sp == "Tat SPI":
            select_analyses = st.radio("Which analyses would you like to take a look at?", ("Full signal peptide length distribution", "Length distributions of each signal peptide region", "Sequence logos"))
            if select_analyses == "Full signal peptide length distribution":
                st.subheader("Histogram of signal peptide length distribution")
                MakeHistograms(tax_id+"_SP_regions_counts.tab", "TAT", group = scient_name, see_regions = False, rank = clean_rank)
            elif select_analyses == "Length distributions of each signal peptide region":
                st.subheader("Histograms of region length distributions")
                MakeHistograms(tax_id+"_SP_regions_counts.tab", "TAT", group = scient_name, see_regions = True, rank = clean_rank)
            elif select_analyses == "Sequence logos":
                st.subheader("Sequence logos aligned at region borders")
                option1,option2 = st.columns(2)
                option1_inf,option2_inf = st.columns(2)
                color_scheme = option1.radio("Select a color scheme for the sequence logos", ("Default", "Hydrophobicity", "Charge", "Chemistry"))
                with option1_inf.expander("Information about color schemes"):
                    st.image(seq_logo_color_schemes)
                    st.write("""
                        Default color scheme divides amino acids into negatively charged (red), 
                        positively charged (blue), hydrophobic (black) and others (green). This color scheme
                        is useful when examining important signal peptides properties such as
                        the positive charge of n-regions and hydrophobicity of h-regions. 
                        """)

                keep_zero = option2.radio("Do you want to include 0 as the position before alignment site?", ("Yes", "No"))
                with option2_inf.expander("Information about position assignments"):
                    st.write("""
                        Choosing either of these options do not influence the logos themselves, 
                        except changing whether the position before alignment site should be 
                        written as position 0 or position -1. 
                        """)
                #st.caption("Alignment at N-H region border")
                min_pos_n = ExtractAlignmentPosition(tax_id, "TAT", "N")[0]
                len_seq_n = ExtractAlignmentPosition(tax_id, "TAT", "N")[1]
                max_pos_n = int(len_seq_n) - int(min_pos_n)
                st.write("Select the range of positions in the logo aligned at the n-h region border")
                low,high = st.columns(2)
                lower_n = low.number_input('Lower value', min_value = -int(min_pos_n), max_value = -1, value = -5, step = 1)
                upper_n = high.number_input('Upper value', min_value = 1, max_value = max_pos_n, value = 7, step = 1)
                MakeLogo(tax_id,"TAT","n", lower_n, upper_n, scient_name, color_scheme, keep_zero)
                #except ValueError:
                #    st.warning("Your selected value for lower part of the range is too small. Please try with a greater value.")
                min_pos_h = ExtractAlignmentPosition(tax_id, "TAT", "H")[0]
                len_seq_h = ExtractAlignmentPosition(tax_id, "TAT", "H")[1]
                max_pos_h = int(len_seq_h) - int(min_pos_h)
                st.write("Select the range of positions in the logo aligned at the h-c region border")
                low,high = st.columns(2)
                lower_h = low.number_input('Lower value', min_value = -int(min_pos_h), max_value = -1, value = -15, step = 1)
                upper_h = high.number_input('Upper value', min_value = 1, max_value = max_pos_h, value = 5, step = 1)
                MakeLogo(tax_id,"TAT","h", lower_h, upper_h, scient_name, color_scheme)

                min_pos_c = ExtractAlignmentPosition(tax_id, "TAT", "C")[0]
                len_seq_c = ExtractAlignmentPosition(tax_id, "TAT", "C")[1]
                max_pos_c = int(len_seq_c) - int(min_pos_c)
                st.write("Select the range of positions in the logo aligned at cleavage site")
                low,high = st.columns(2)
                lower_c = low.number_input('Lower value', min_value = -int(min_pos_c), max_value = -1, value = -7, step = 1)
                upper_c = high.number_input('Upper value', min_value = 1, max_value = max_pos_c, value = 5, step = 1)
                MakeLogo(tax_id,"TAT","c", lower_c, upper_c, scient_name, color_scheme, keep_zero)

                with st.expander("Creation of the sequence logos"):
                    st.write("""
                        Write something about how the logos are created (Shannon, Schneider and Steven, 1990).
                        """)
                with st.expander("Calculation of errorbars"):
                    st.write("""
                        The magnitude of the errorbars are based on the small sample correction
                        factor. The small sample correction factor for sequence logos with 
                        amino acids is calculated as
                        """)
                    st.latex(r'''e_n = \frac{1}{\ln{2}} \cdot \frac{20-1}{2\cdot n}''')
                    st.write("""
                        The height of the errorbars correspond to twice the small sample correction
                        factor. 
                        """)
        
        #Tat SPII analyses
        elif select_sp == "Tat SPII":
            select_analyses = st.radio("Which analyses would you like to take a look at?", ("Full signal peptide length distribution", "Length distributions of each signal peptide region", "Sequence logos"))
            if select_analyses == "Full signal peptide length distribution":
                st.subheader("Histogram of signal peptide length distribution")
                MakeHistograms(tax_id+"_SP_regions_counts.tab", "TATLIPO", group = scient_name, see_regions = False, rank = clean_rank)
            elif select_analyses == "Length distributions of each signal peptide region":
                st.subheader("Histograms of region length distributions")
                MakeHistograms(tax_id+"_SP_regions_counts.tab", "TATLIPO", group = scient_name, see_regions = True, rank = clean_rank)
            elif select_analyses == "Sequence logos":
                st.subheader("Sequence logos aligned at region borders")
                option1,option2 = st.columns(2)
                option1_inf,option2_inf = st.columns(2)
                color_scheme = option1.radio("Select a color scheme for the sequence logos", ("Default", "Hydrophobicity", "Charge", "Chemistry"))
                with option1_inf.expander("Information about color schemes"):
                    st.image(seq_logo_color_schemes)
                    st.write("""
                        Default color scheme divides amino acids into negatively charged (red), 
                        positively charged (blue), hydrophobic (black) and others (green). This color scheme
                        is useful when examining important signal peptides properties such as
                        the positive charge of n-regions and hydrophobicity of h-regions. 
                        """)

                keep_zero = option2.radio("Do you want to include 0 as the position before alignment site?", ("Yes", "No"))
                with option2_inf.expander("Information about position assignments"):
                    st.write("""
                        Choosing either of these options do not influence the logos themselves, 
                        except changing whether the position before alignment site should be 
                        written as position 0 or position -1. 
                        """)

                #st.caption("Alignment at N-H region border")
                min_pos_n = ExtractAlignmentPosition(tax_id, "TATLIPO", "N")[0]
                len_seq_n = ExtractAlignmentPosition(tax_id, "TATLIPO", "N")[1]
                max_pos_n = int(len_seq_n) - int(min_pos_n)
                st.write("Select the range of positions in the logo aligned at the n-h region border")
                low,high = st.columns(2)
                lower_n = low.number_input('Lower value', min_value = -int(min_pos_n), max_value = -1, value = -5, step = 1)
                upper_n = high.number_input('Upper value', min_value = 1, max_value = max_pos_n, value = 7, step = 1)
                MakeLogo(tax_id,"TATLIPO","n", lower_n, upper_n, scient_name, color_scheme, keep_zero)
                #except ValueError:
                #    st.warning("Your selected value for lower part of the range is too small. Please try with a greater value.")
                min_pos_h = ExtractAlignmentPosition(tax_id, "LIPO", "H")[0]
                len_seq_h = ExtractAlignmentPosition(tax_id, "LIPO", "H")[1]
                max_pos_h = int(len_seq_h) - int(min_pos_h)
                st.write("Select the range of positions in the logo aligned at the h-c region border")
                low,high = st.columns(2)
                lower_h = low.number_input('Lower value', min_value = -int(min_pos_h), max_value = -1, value = -7, step = 1)
                upper_h = high.number_input('Upper value', min_value = 1, max_value = max_pos_h, value = 5, step = 1)
                MakeLogo(tax_id,"TATLIPO","h", lower_h, upper_h, scient_name, color_scheme, keep_zero)

                with st.expander("Creation of the sequence logos"):
                    st.write("""
                        Write something about how the logos are created (Shannon, Schneider and Steven, 1990).
                        """)
                with st.expander("Calculation of errorbars"):
                    st.write("""
                        The magnitude of the errorbars are based on the small sample correction
                        factor. The small sample correction factor for sequence logos with 
                        amino acids is calculated as
                        """)
                    st.latex(r'''e_n = \frac{1}{\ln{2}} \cdot \frac{20-1}{2\cdot n}''')
                    st.write("""
                        The height of the errorbars correspond to twice the small sample correction
                        factor. 
                        """)


    except NameError as err:
        st.warning('The program was not able to evoke analyses of ' + option + '. Please try with another taxonomical group or organism name.')
        st.write(str(err))

