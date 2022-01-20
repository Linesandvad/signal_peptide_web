#!/usr/bin/env python3

import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import math
import streamlit as st
import sys
from PIL import Image
import logomaker

def ExtractAlignmentPosition(tax_id,sp_type,region):
	os.chdir("./Aligned_sequences/"+tax_id)

	region_lowercase = region.lower()
	region_uppercase = region.upper()

	sp_type = sp_type.upper()

	infile = open(sp_type+"_align_"+region_lowercase+".txt", "r")
	line = infile.readline()
	len_seqs = len(line[:-1])

	infile.close()

	infile_positions = open("align_positions.tab", "r")
	region_matched = False

	try:

		for line in infile_positions:
			splits = line.split("\t")

			if splits[0] == sp_type:
				if splits[1] == region_uppercase:
					alignment_position = splits[2] 		#Set position 1 of logo (selected region border)
					region_matched = True

		if not region_matched:
			raise ValueError
	except ValueError as err:
		print("Error! The signal peptide type {0} does not have an {1}-region.".format(sp_type,region))

	infile_positions.close()

	os.chdir("../../")

	return alignment_position, len_seqs

def MakeLogo(tax_id,sp_type,region, low_limit, high_limit, scient_name, color_scheme, keep_zero):
	#Define color scheme

	color_1990 = {"K": "blue", "R": "blue", "H": "blue",
	"D": "red", "E": "red",
	"A": "black", "V": "black", "L": "black", "I": "black",
	"P": "black", "W": "black", "F": "black", "M": "black",
	"N": "green", "G": "green", "Q": "green", "S": "green",
	"C": "green", "Y": "green", "T": "green",
	"X": "magenta"}

	if color_scheme == "Default":
		color_logo = color_1990
	elif color_scheme == "Hydrophobicity":
		color_logo = "hydrophobicity"
	elif color_scheme == "Charge":
		color_logo = "charge"
	elif color_scheme == "Chemistry":
		color_logo = "chemistry"

	region_lowercase = region.lower()
	region_uppercase = region.upper()

	sp_type = sp_type.upper()

	os.chdir("./Aligned_sequences/"+tax_id)
	infile = open(sp_type+"_align_"+region_lowercase+".txt", "r")

	seqs = []

	#Store all SP sequences as strings in a list
	for line in infile:
		seqs.append(line[:-1])

	infile.close()


	infile_positions = open("align_positions.tab", "r")
	region_matched = False

	try:

		for line in infile_positions:
			splits = line.split("\t")

			if splits[0] == sp_type:
				if splits[1] == region_uppercase:
					alignment_position = splits[2] 		#Set position 1 of logo (selected region border)
					region_matched = True

		if not region_matched:
			raise ValueError
	except ValueError as err:
		print("Error! The signal peptide type {0} does not have an {1}-region.".format(sp_type,region))

	infile_positions.close()


	if int(alignment_position) == 0:
		st.info("No signal peptides of type " + sp_type + " were detected for this group.")
		os.chdir("../../")

	else:

		def __SetPositions(low_limit, high_limit, keep_zero):
			try:
				if int(low_limit) > 0:
					raise ValueError
			except ValueError as err:
				print("The lower limit of the selected range must be a negative number. Your input lower input number {0} will be converted to -{0}.".format(low_limit))
				low_limit = -int(low_limit)
			try:
				if int(high_limit) < 0:
					raise ValueError
			except ValueError as err:
				print("The high limit of the selected range must be a positive number. Your input upper input number {0} will be converted to {1}.".format(str(high_limit), str(high_limit)[1:]))
				high_limit = str(high_limit)[1:]
				high_limit = int(high_limit)

			positions = list(range(low_limit,high_limit+1))

			if keep_zero == "Yes":
				positions.remove(int(low_limit))
			elif keep_zero == "No":
				positions.remove(0)

			number_of_pos = int(high_limit) - int(low_limit)

			set_boundary = -low_limit - 0.5

			return positions, low_limit, high_limit, number_of_pos, set_boundary

		func_output = __SetPositions(low_limit, high_limit, keep_zero)
		pos_range = func_output[0]
		lower = func_output[1] + int(alignment_position)
		upper = func_output[2] + int(alignment_position) - 1
		number_of_pos = func_output[3]
		set_boundary = func_output[4]


		#Compute small smaple corrections
		count_df = logomaker.alignment_to_matrix(sequences = seqs, to_type = "counts", characters_to_ignore = "-")

		Rseq = []
		inf_df = count_df.copy()

		for i in range(len(count_df)):
			n = sum(count_df.loc[i])
			if n != 0:
				e_n = (20-1)/(2*math.log(2)*n)
				Hl = 0
				for j in range(len(count_df.loc[i])):
					freq = (count_df.loc[i][j]/n)
					if freq > 0:
						logfreq = math.log2(freq)
					else:
						logfreq = 0
					Hl += freq*logfreq
			else:
				e_n = 0
				Hl = 0
				freq = 0
				logfreq = 0
			H = -Hl
			if math.log2(20)-(H+e_n) < 0:
				Rseq.append(0)
			else:
				Rseq.append(math.log2(20)-(H+e_n))

		for i in range(len(count_df)):
			n = sum(count_df.loc[i])
			for j in range(len(count_df.loc[i])):
				if n != 0:
					freq = count_df.loc[i][j]/n
					inf_df.loc[i][j] = freq*Rseq[i]
				else:
					freq = 0
					inf_df.loc[i][j] = freq*Rseq[i]

		inf_df = inf_df.loc[lower:upper]
		inf_df_cleaned = inf_df.copy()

		inf_df_cleaned.index = list(range(0,number_of_pos))


		count_df = count_df.loc[lower:upper]
		count_df_cleaned = count_df.copy()
		count_df_cleaned.index = list(range(0,number_of_pos))

		E_n = []
		heights = []

		for i in range(len(count_df_cleaned)):
		    n = sum(count_df_cleaned.loc[i])
		    if n != 0:
		    	e_n = (20-1)/(2*math.log(2)*n)
		    
		    else:
		    	e_n = 0
		    errorbar = 2*e_n
		    E_n.append(errorbar)

		    heights.append(sum(inf_df_cleaned.loc[i]))



		##Preparation for logo styling
		if sp_type == "SP":
			sp_type_desc = "Sec SPI"
			if region_lowercase == "n":
				region_border = "n-h region border"
			elif region_lowercase == "h":
				region_border = "h-c region border"
			elif region_lowercase == "c":
				region_border = "cleavage site"
		
		elif sp_type == "LIPO":
			sp_type_desc = "Sec SPII"
			if region_lowercase == "n":
				region_border = "n-h region border"
			elif region_lowercase == "h":
				region_border = "cleavage site"

		elif sp_type == "PILIN":
			sp_type_desc = "Sec SPIII"
			if region_lowercase == "n":
				region_border = "cleavage site"
		
		elif sp_type == "TAT":
			sp_type_desc = "Tat SPI"
			if region_lowercase == "n":
				region_border = "n-h region border"
			elif region_lowercase == "h":
				region_border = "h-c region border"
			elif region_lowercase == "c":
				region_border = "cleavage site"
		
		elif sp_type == "TATLIPO":
			sp_type_desc = "Tat SPII"
			if region_lowercase == "n":
				region_border = "n-h region border"
			elif region_lowercase == "h":
				region_border = "cleavage site"

		##Make the Logo
		#fig, ax = plt.subplots(1,1,figsize=[4,2])
		#fig = plt.figure(figsize=[6, 2])

		fig, ax = plt.subplots(1,1,figsize=[4,2])

		#seq_logo = plt.figure(figsize=[8, 2])
		seq_logo = logomaker.Logo(inf_df_cleaned,
		stack_order = "big_on_top",
		font_name='DejaVu Sans',
		color_scheme=color_logo,
		vpad=0)

		# style using Logo methods
		seq_logo.style_xticks(anchor=0, spacing=30, rotation=0)

		seq_logo.ax.set_xticks(range(len(pos_range)))
		seq_logo.ax.set_xticklabels(x for x in pos_range)
		seq_logo.ax.axvline(set_boundary, color='k', linewidth=1, linestyle=':')
		seq_logo.ax.set_title(sp_type_desc + " " + region_border + " in " + scient_name)
		seq_logo.ax.set_ylabel('Information (bits)')
		seq_logo.ax.set_xlabel('Position')
		seq_logo.ax.set_xlim([-1, len(inf_df_cleaned)])
		seq_logo.ax.set_ylim([0, math.log2(20)])


		x = np.linspace(0,number_of_pos-1,number_of_pos)
		y = heights

		xerr = 0
		yerr = E_n
		 
		seq_logo.ax.errorbar(x, y,
		            xerr=xerr,
		            yerr=yerr,
		            fmt='.',
		           color='grey',
		           ecolor='grey')
		seq_logo.ax.grid(False)


		os.chdir("../../Logos_saved")

		plt.savefig("seq_logo.png")
		image = Image.open('seq_logo.png')

		os.chdir("../")

		st.image(image)

		#plt.show()


#Example shows
#MakeLogo("1783275","sp","n", -5, 10, "test")

##Handle this:
#MakeLogo("2158","tatlipo","n", -5, 10)

#And title
