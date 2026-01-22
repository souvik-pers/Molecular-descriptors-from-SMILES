#@author: Souvik Chakraborty, IHPC, Singapore

#https://github.com/ecrl/padelpy

from padelpy import from_smiles
import numpy as np
import pandas as pd


#"""
#First import data from the text file using NumPy's genfromtxt
smiles = np.genfromtxt("./New_smiles.csv",delimiter="\t", dtype=None, comments=None, encoding=None, names=None)
#print(smiles)
#print(len(smiles))


smiles_list = [x for x in smiles.tolist()]
print(smiles_list)
print(len(smiles_list))


##SERIAL INPUT
i=0
for item in smiles_list:
    #print(item)
    i=i+1
    descriptors = from_smiles(item, output_csv='./descriptors-out-%s.csv' %i, threads = 10, timeout=360)
    #Appending SMILES entry at the LAST column
    df=pd.read_csv('./descriptors-out-%s.csv' %i)
    df["smiles"]=item
    df.to_csv('./descriptors-out-%s.csv' %i)
    

## Next Step: To combine all individual output csv files in a single file

import glob
import os

path = r'./' # use your path
all_files = glob.glob(os.path.join(path, "descriptors-out-*.csv"))
#print(all_files)
#print(len(all_files))

#Need sorting filename's order in all_files list
#https://copyprogramming.com/howto/python-sort-file-names-with-numbers
import natsort
all_files_sorted = natsort.natsorted(all_files)
print(all_files_sorted)


df = pd.concat((pd.read_csv(f) for f in all_files_sorted), ignore_index=True)
print(df.shape)
df.to_csv('./New_smiles_PADEL_FEATURES.csv')
