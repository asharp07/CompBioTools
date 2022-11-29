# Written by Amanda K. Sharp
# In the Brown Lab
#
# Version 2
# Last updated 3/18/21

import pandas as pd 
import numpy as np
import os
import Interactions as inter

def gather_finger_data(residue_key, my_dir, interaction, datatype, ignorechain):
    finger_data = []
    ligands = []
    files = [f for f in os.listdir('.') if os.path.isfile(f)]
    for file in files:
        if '.csv' in file and 'fingerprint' in file:
            protein, ligand, ext = file.split('_')
            ligands.append(ligand)
            pd.read_csv(file, header=None).T.to_csv(f'temp_{protein}_{ligand}_flip.csv', header=False, index=False)
            df = pd.read_csv(f'temp_{protein}_{ligand}_flip.csv')
            col_names = []
            num_poses = []
            for i in range(len(df.columns)):
                if i == 0:
                    col_names.append('Title')
                else:
                    col_names.append(f'Pose{i}')
                    num_poses.append(i)
            df.columns = col_names
            new = df['Title'].str.split('_', n=1, expand=True)
            df['Chain_Res'] = new[0]
            df['interaction type'] = new[1]
            df['Residue_Number'] = df['Chain_Res'].str.replace('([A-Z]+)', '')
            df.drop(df.tail(3).index,inplace=True)
            df_updated = append_residue(df, residue_key, protein, ligand, num_poses, interaction, datatype, ignorechain)
            finger_data.append(df_updated.values)
            os.remove(f'temp_{protein}_{ligand}_flip.csv')
    return finger_data, ligands

def append_residue(df, residue_key, protein, ligand, num_poses, interaction, datatype, ignorechain):
    int_dict = {'all' : inter.All_Contact, 
        'charged' : inter.Charged, 
        'polar' : inter.Polar, 
        'hydrophobic' : inter.Hydrophobic,
        'backbone' : inter.Backbone}
    try:
        data = int_dict[interaction](df, num_poses)
    except KeyError:
        print('Error in Interaction Type')
    data = data.reset_index()
    residue_key = residue_key.reset_index()
    res_match = pd.concat([data, residue_key], ignore_index=True, axis='columns')
    pose_col = ['index1', 'Title']
    for i in num_poses:
        pose_col.append(f'Pose{i}')
    end_col = ['Chain_Res', 'interaction type', 'Res',  'total',  'normalized', 'index', 'Res_Chain_Number', 'Residue_Name']
    res_match.columns = pose_col + end_col
    res_match['Residue_Name'] = res_match['Residue_Name'].str.replace('^ +| +$', '')
    res_match['ChainID'] = res_match['Res_Chain_Number'].str.extract('([A-Za-z]+)', expand = True)
    res_match['Res_Number'] = res_match['Res_Chain_Number'].str.replace('([A-Za-z]+)', '')
    if ignorechain == 'yes':
        res_match['label'] = res_match[['Res_Number', 'Residue_Name']].apply(lambda row: ''.join(row.values.astype(str)), axis=1)
    else:
        res_match['label_pt1'] = res_match[['Res_Number', 'ChainID']].apply(lambda row: ''.join(row.values.astype(str)), axis=1)
        res_match['label'] = res_match[['label_pt1', 'Residue_Name']].apply(lambda row: '-'.join(row.values.astype(str)), axis=1)
    data = res_match[[f'label', f'{datatype}']]
    return data.transpose()
