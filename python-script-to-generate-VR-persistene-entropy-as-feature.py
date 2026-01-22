#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 15:20:28 2023

@author: souvik
"""

# RIPSER link
# https://github.com/Ripser/ripser
# Ripser is a lean C++ code for the computation of Vietoris–Rips persistence barcodes.
# It can do just this one thing, but does it extremely well.

# https://ripser.scikit-tda.org/en/latest/reference/stubs/ripser.ripser.html

import pandas as pd
import os, os.path
#import natsort
from sklearn.metrics import pairwise_distances
from ripser import ripser
from persim import plot_diagrams
import matplotlib.pyplot as plt
import pymatgen.core as pmg
#import ripser
from persim.persistent_entropy import *

#################  PERSISTENCE ENTROPY  #######################

#https://persim.scikit-tda.org/en/latest/notebooks/Persistence%20barcode%20measure.html

inputfolder = "./split_xyz"  # path for dir
lst = os.listdir(inputfolder) # your directory path
num_files = len(lst)
print(num_files)


i=0
for index in range(1,num_files+1):
    i=i+1
    incr_filename = inputfolder + '/molecule_' + str(index) + '.xyz'
    print(incr_filename)
    
    molecule = pmg.Molecule.from_file(incr_filename)
    # Create a distance matrix based on neighbor information
    distances = pairwise_distances(molecule.cart_coords, metric="euclidean")

    # Generate the Vietoris-Rips complex using Ripser
    result = ripser(distances, maxdim=2)  # You can change maxdim to adjust the complex dimension

    #print(result)
    
    entropy_list_in_all_dim = []
    for j in range(3):
        diag = result['dgms'][j]
        entropy = persistent_entropy(diag)
        entropy_list_in_all_dim.extend(entropy)
        
    print('entropy_in_three_dimensions = ', entropy_list_in_all_dim)
    print()
    print()
    
    
    
    df = pd.DataFrame([list(entropy_list_in_all_dim)])
    #print(df)
    # Define the number of columns and a prefix for header names
    num_columns = len(list(entropy_list_in_all_dim))
    header_prefix = "PERSISTENCE_ENTROPY_DIM_"

    # Create a list to store the header names
    headers = [f"{header_prefix}{i}" for i in range(0, num_columns)]

    # Print the list of header names
    #print(headers)

    # Assign the headers to the DataFrame
    #df.columns = headers

    # Print the DataFrame
    #print(df)
    # Assign the headers to the DataFrame
    df.columns = headers

    # Print the DataFrame
    #print(df)
    df.to_csv('./persistence-entropy-out-%s.csv' %i)
    
    
print()

## Next Step: To combine all individual output csv files in a single file
#https://stackoverflow.com/questions/20906474/import-multiple-csv-files-into-pandas-and-concatenate-into-one-dataframe

import glob
import os


path = r'./' # use your path
all_files = glob.glob(os.path.join(path, "./persistence-entropy-out-*.csv"))
#print(all_files)
#print(len(all_files))

#Need sorting filename's order in all_files list
#https://copyprogramming.com/howto/python-sort-file-names-with-numbers
import natsort
all_files_sorted = natsort.natsorted(all_files)
print(all_files_sorted)


df = pd.concat((pd.read_csv(f) for f in all_files_sorted), ignore_index=True)
print(df.shape)
df.to_csv('./New_smiles_VIETORIS_RIPS_PERSISTNCE_ENTROPY_OUT.csv')
