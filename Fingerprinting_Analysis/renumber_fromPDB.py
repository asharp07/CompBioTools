# Written by Amanda K. Sharp
# In the Brown Lab
#
# Version 2
# Last updated 3/18/21

import csv
import pandas as pd 
import os

def Protein_Key(Protein, my_dir):
    with open(Protein, "rt", newline='') as input_file:
        res_info = []
        chain_resnum = []
        for line in input_file:
            if 'ATOM' in line:
                chain_resnum.append(line[21:26])
                res_info.append(line[17:20])
        protein = dict(zip(chain_resnum, res_info))
    df = pd.DataFrame(list(protein.items()))
    df.columns = ['Res_Chain_Number', 'Residue_Name']
    return df
