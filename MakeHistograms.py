#!/usr/bin/env python3

import streamlit as st
import matplotlib.pyplot as plt
import math
from CountSPs import *


def MakeHistograms(filename, sp_type, group, see_regions, rank):
    """Make histograms of SP region length distributions, given a phylogenetic group and SP type"""

    #Collect data from CountSPs
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

    #Collect full SP length data
    len_sp = []
    len_lipo = []
    len_tat = []
    len_tatlipo = []
    len_pilin = []

    #Sec SPI
    for i in range(len(n_region_sp)):
        sp_len = int(n_region_sp[i]) + int(h_region_sp[i]) + int(c_region_sp[i])
        len_sp.append(sp_len)

    #Sec SPII
    for i in range(len(n_region_lipo)):
        lipo_len = int(n_region_lipo[i]) + int(h_region_lipo[i])
        len_lipo.append(lipo_len)

    #Tat SPI
    for i in range(len(n_region_tat)):
        tat_len = int(n_region_tat[i]) + int(h_region_tat[i]) + int(c_region_tat[i])
        len_tat.append(tat_len)

    #Tat SPII
    for i in range(len(n_region_tatlipo)):
        tatlipo_len = int(n_region_tatlipo[i]) + int(h_region_tatlipo[i])
        len_tatlipo.append(tatlipo_len)

    #Sec SPIII
    for i in range(len(n_region_pilin)):
        pilin_len = int(n_region_pilin[i])
        len_pilin.append(pilin_len)

    
    def __Histograms(region, sp_type, group, region_lengths, w, hist_col, rank):
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
            ax.set_title("Distribution of "+ sp_type + " " + region + "-region lengths \n for "+ rank + " " + group)
            #Display histogram in streamlit app
            st.pyplot(fig)

    #Use this function if all region counts in a region for a group is the same (rare case)
    def __HistogramsOneLength(region, sp_type, group, region_lengths, hist_col, rank):
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
            ax.set_title("Distribution of "+ sp_type + " " + region + "-region lengths \n for "+ rank + " " + group)
            #Display histogram in streamlit app
            st.pyplot(fig)

    def __HistogramsCompleteSP(sp_type, group, sp_lengths, w, hist_col, rank):
        """Make histograms, display to streamlit"""
        #Sort region lengths from shortest to longest
        sp_lengths.sort()
        #Only produce non-empty histograms
        if len(sp_lengths) != 0: 
            #Define number of bins (1 bin for each integer)
            n = math.ceil((sp_lengths[-1] - sp_lengths[0])/w)
            plt.style.use('ggplot')
            fig, ax = plt.subplots()
            #Create histogram
            ax.hist(sp_lengths, bins = n, density = True, color = hist_col)
            #Create histogram title
            ax.set_title("Distribution of complete "+ sp_type + " lengths \n for "+ rank + " " + group)
            #Display histogram in streamlit app
            st.pyplot(fig)

    def __HistogramsCompleteSPOneLength(sp_type, group, sp_lengths, hist_col, rank):
        """Make histograms, display to streamlit"""
        #Sort region lengths from shortest to longest
        sp_lengths.sort()
        #Only produce non-empty histograms
        if len(sp_lengths) != 0: 
            #Define number of bins (1 bin for each integer)
            n = 1
            plt.style.use('ggplot')
            fig, ax = plt.subplots()
            #Create histogram
            ax.hist(sp_lengths, bins = n, density = True, color = hist_col)
            #Create histogram title
            ax.set_title("Distribution of complete "+ sp_type + " lengths \n for "+ rank + " " + group)
            #Display histogram in streamlit app
            st.pyplot(fig)

    #Display histograms of SP region lengths of the selected SP type
    #Displaying histograms for Sec SPI
    if see_regions:
        if sp_type == "SP":
            if len(n_region_sp) == 0:
                st.info("No signal peptides of the chosen type were found for this group.")
            else:
                try:
                    __Histograms("n", "sec SPI", group, n_region_sp, 1, "blue", rank)
                except ValueError as err:
                    __HistogramsOneLength("n", "sec SPI", group, n_region_sp, "blue", rank)
                    if len(n_region_sp) > 1:
                        st.info("All Sec SPI n-regions within this group consist of " + str(n_region_sp[0]) + " residues.")
                try:
                    __Histograms("h", "sec SPI", group, h_region_sp, 1, "green", rank)
                except ValueError as err:
                    __HistogramsOneLength("h", "sec SPI", group, h_region_sp, "green", rank)
                    if len(h_region_sp) > 1:
                        st.info("All Sec SPI h-regions within this group consist of " + str(h_region_sp[0]) + " residues.")
                try:
                    __Histograms("c", "sec SPI", group, c_region_sp, 1, "red", rank)
                except ValueError as err:
                    __HistogramsOneLength("c", "sec SPI", group, c_region_sp, "red", rank)
                    if len(c_region_sp) > 1:
                        st.info("All Sec SPI c-regions within this group consist of " + str(c_region_sp[0]) + " residues.")
        
        #Displaying histograms for Sec SPII
        elif sp_type == "LIPO":
            if len(n_region_lipo) == 0:
                st.info("No signal peptides of the chosen type were found for this group.")
            else:
                try:
                    __Histograms("n", "sec SPII", group, n_region_lipo, 1, "blue", rank)
                except ValueError as err:
                    __HistogramsOneLength("n", "sec SPII", group, n_region_lipo, "blue", rank)
                    if len(n_region_lipo) > 1:
                        st.info("All Sec SPII n-regions within this group consist of " + str(n_region_lipo[0]) + " residues.")

                try:
                    __Histograms("h", "sec SPII", group, h_region_lipo, 1, "green", rank)
                except ValueError as err:
                    __HistogramsOneLength("h", "sec SPII", group, h_region_lipo, "green", rank)
                    if len(h_region_lipo) > 1:
                        st.info("All Sec SPII h-regions within this group consist of " + str(h_region_lipo[0]) + " residues.")
        
        #Displaying histograms for Tat SPI
        elif sp_type == "TAT":
            if len(n_region_tat) == 0:
                st.info("No signal peptides of the chosen type were found for this group.")
            else:
                try:
                    __Histograms("n", "tat SPI", group, n_region_tat, 1, "blue", rank)
                except ValueError as err:
                    __HistogramsOneLength("n", "tat SPI", group, n_region_tat, "blue", rank)
                    if len(n_region_tat) > 1:
                        st.info("All Tat SPI n-regions within this group consist of " + str(n_region_tat[0]) + " residues.")
                try:
                    __Histograms("h", "tat SPI", group, h_region_tat, 1, "green", rank)
                except ValueError as err:
                    __HistogramsOneLength("h", "tat SPI", group, h_region_tat, "green", rank)
                    if len(h_region_tat) > 1:
                        st.info("All Tat SPI h-regions within this group consist of " + str(h_region_tat[0]) + " residues.")
                try:
                    __Histograms("c", "tat SPI", group, c_region_tat, 1, "red", rank)
                except ValueError as err:
                    __HistogramsOneLength("c", "tat SPI", group, c_region_tat, "red", rank)
                    if len(c_region_tat) > 1:
                        st.info("All Tat SPI c-regions within this group consist of " + str(c_region_tat[0]) + " residues.")
        
        #Displaying histograms for Tat SPII
        elif sp_type == "TATLIPO":
            if len(n_region_tatlipo) == 0:
                st.info("No signal peptides of the chosen type were found for this group.")
            else:
                try:
                    __Histograms("n", "tat SPII", group, n_region_tatlipo, 1, "blue", rank)
                except ValueError as err:
                    __HistogramsOneLength("n", "tat SPII", group, n_region_tatlipo, "blue", rank)
                    if len(n_region_tatlipo) > 1:
                        st.info("All Tat SPII n-regions within this group consist of " + str(n_region_tatlipo[0]) + " residues.")
                try:
                    __Histograms("h", "tat SPII", group, h_region_tatlipo, 1, "green", rank)
                except ValueError as err:
                    __HistogramsOneLength("h", "tat SPII", group, h_region_tatlipo, "green", rank)
                    if len(h_region_tatlipo) > 1:
                        st.info("All Tat SPII h-regions within this group consist of " + str(h_region_tatlipo[0]) + " residues.")
        
        #Displaying histograms for Sec SPIII
        elif sp_type == "PILIN":
            if len(n_region_pilin) == 0:
                st.info("No signal peptides of the chosen type were found for this group.")
            else:
                try:
                    __Histograms("n", "sec SPIII", group, n_region_pilin, 1, "blue", rank)
                except ValueError as err:
                    __HistogramsOneLength("n", "sec SPIII", group, n_region_pilin, "blue", rank)
                    if len(n_region_pilin) > 1:
                        st.info("All Sec SPIII signal peptides within this group contains " + str(n_region_pilin[0]) + " residues.")

    if not see_regions:
        if sp_type == "SP":
            if len(len_sp) == 0:
                st.info("No signal peptides of the chosen type were found for this group.")
            else:
                try:
                    __HistogramsCompleteSP("sec SPI", group, len_sp, 1, "purple", rank)
                except ValueError as err:
                    __HistogramsCompleteSPOneLength("sec SPI", group, len_sp, "purple", rank)
                    if len(len_sp) > 1:
                        st.info("All Sec SPI signal peptides within this group consist of " + str(len_sp[0]) + " residues.")
        #Displaying histograms for Sec SPII
        elif sp_type == "LIPO":
            if len(len_lipo) == 0:
                st.info("No signal peptides of the chosen type were found for this group.")
            else:
                try:
                    __HistogramsCompleteSP("sec SPII", group, len_lipo, 1, "purple", rank)
                except ValueError as err:
                    __HistogramsCompleteSPOneLength("sec SPII", group, len_lipo, "purple", rank)
                    if len(len_lipo) > 1:
                        st.info("All Sec SPII signal peptides within this group consist of " + str(len_lipo[0]) + " residues.")

        #Displaying histograms for Tat SPI
        elif sp_type == "TAT":
            if len(len_tat) == 0:
                st.info("No signal peptides of the chosen type were found for this group.")
            else:
                try:
                    __HistogramsCompleteSP("tat SPI", group, len_tat, 1, "purple", rank)
                except ValueError as err:
                    __HistogramsCompleteSPOneLength("tat SPI", group, len_tat, "purple", rank)
                    if len(len_tat) > 1:
                        st.info("All Tat SPI signal peptides within this group consist of " + str(len_tat[0]) + " residues.")

        #Displaying histograms for Tat SPII
        elif sp_type == "TATLIPO":
            if len(len_tatlipo) == 0:
                st.info("No signal peptides of the chosen type were found for this group.")
            else:
                try:
                    __HistogramsCompleteSP("tat SPII", group, len_tatlipo, 1, "purple", rank)
                except ValueError as err:
                    __HistogramsCompleteSPOneLength("tat SPII", group, len_tatlipo, "purple", rank)
                    if len(len_tatlipo) > 1:
                        st.info("All Tat SPII signal peptides within this group consist of " + str(len_tatlipo[0]) + " residues.")


        #Displaying histograms for Sec SPIII
        elif sp_type == "PILIN":
            if len(len_pilin) == 0:
                st.info("No signal peptides of the chosen type were found for this group.")
            else:
                try:
                    __HistogramsCompleteSP("sec SPIII", group, len_pilin, 1, "purple", rank)
                except ValueError as err:
                    __HistogramsCompleteSPOneLength("sec SPIII", group, len_pilin, "purple", rank)
                    if len(len_pilin) > 1:
                        st.info("All Sec SPIII signal peptides within this group consist of " + str(len_pilin[0]) + " residues.")
