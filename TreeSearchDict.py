#!/usr/bin/env python3

import pickle

def TreeSearchDict(domain):
    infile_dict = open(domain+"_tax_groups_children.pk", "rb")
    tree_search_file = pickle.load(infile_dict)

    return tree_search_file
