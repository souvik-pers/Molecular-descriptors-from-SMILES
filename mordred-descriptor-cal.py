#Coded by Souvik Chakraborty

from rdkit import Chem
from mordred import Calculator, descriptors
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

#https://github.com/mordred-descriptor/mordred
calc = Calculator(descriptors, ignore_3D=True)

mols = [Chem.MolFromSmiles(smi) for smi in smiles_list]

df = calc.pandas(mols)

print (df.shape)

df.to_csv('./New_smiles-Mordred-descriptor-out.csv')

