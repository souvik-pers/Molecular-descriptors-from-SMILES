#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: souvik chakraborty, IHPC, Singapore
"""

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
import numpy as np



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
    rdkit_2d_desc = []
    i=i+1
    calc = MoleculeDescriptors.MolecularDescriptorCalculator([x[0] for x in Descriptors._descList])
    header = calc.GetDescriptorNames()
    mol = Chem.MolFromSmiles(item)
    
    
    
    ds = calc.CalcDescriptors(mol)
    rdkit_2d_desc.append(ds)
    df = pd.DataFrame(rdkit_2d_desc, columns=header)
    df.to_csv('./descriptors-out-%s.csv' %i)
    #Appending SMILES entry at the LAST column
    df=pd.read_csv('./descriptors-out-%s.csv' %i)
    df["smiles"]=item
    df.to_csv('./descriptors-out-%s.csv' %i)
    
    
print()   
#"""


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
df.to_csv('./New_smiles_RDKIT_FEATURES.csv')