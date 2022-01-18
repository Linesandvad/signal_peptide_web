#!/usr/bin/env python3

import streamlit as st
import sys
import os
from ExtractProtCounts import *


def VisualizeProtCounts(tax_id):
    """Display protein- and proteome counts for a given taxonomical group"""
    #Save protein and proteome counts
    prot_count = ExtractProtCounts(tax_id+"_SP_types_counts.tab")[0]
    proteome_count = ExtractProtCounts(tax_id+"_SP_types_counts.tab")[1]

    proteome, protein = st.columns(2)

    proteome.metric("Proteome count", proteome_count) 
    protein.metric("Protein count", prot_count) 