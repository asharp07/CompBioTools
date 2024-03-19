# Written by Amanda K. Sharp
# In the Brown Lab
#
# Version 1
# Last updated 8/7/2023


import numpy as np
import os
import pandas as pd
import itertools
import matplotlib.pyplot as plt
import matplotlib
from scipy.spatial.distance import pdist
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from collections import Counter
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit import RDLogger  
RDLogger.DisableLog('rdApp.*')


class Protein_Analysis:

    # This class processes pdb input for fingerprinting analysis. This also has built in features to export information 
    # about the protein structure
    
    def __init__(self, prot, ligand):
        self.ligand_mol = ligand
        self.prot = Initialize_PDB(prot)
        self.ligand = Initialize_PDB(ligand)
        
    
    def find_atom(res, atom_list, atom_sel=None, ligand=False):

        # Finds individual atom XYZ coordinates 

        atoms = []
        if atom_sel != None:
            for atom in atom_list:
                for sel in atom_sel:
                    if sel in atom:
                        atoms.append(atom)
        
        else:
            atoms = atom_list

        res = res.reset_index()
        atom_coors = {}

        for atom_type in atoms:
            if ligand == True:
                mut = res[res['ATOM_NAME_NUM'] == atom_type]
            else:
                mut = res[res['ATOM_TYPE'] == atom_type]
            atom_coors[atom_type] = mut[['X', 'Y', 'Z']].values.tolist()[0]

        return atom_coors
    
    def split_poses(file):
        poses = []
        temp_data = []

        with open(file, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if line.startswith('ENDMDL'):
                    if len(temp_data) > 0:
                        poses.append(temp_data)
                        temp_data=[]
                elif 'HETATM' in line or 'ATOM' in line:
                    temp_data.append(line)

        output = file.split('.')[0]
        
        if not os.path.exists(output):
            os.makedirs(output)

        for i in range(1, len(poses)+1):
            new_file = open(f'{output}/{output}_{i}.pdb', 'w+')
            for line in poses[i-1]:
                new_file.write(line)
    
    def fingerprinting(self, proximity=False, include_nones=False, atom_num=False, atoms=None, reslabel='multiple', cutoff=5, smiles=None):

        # This function preprocess the pdb file to be organized then anayzed for protein-ligand interactions
        # Options - proximity allows for atoms within close proximity but not an interaction to be listed. These are often conflicts.
        # Inlcude nones is useful when comparing multiple ligands interacting with the same protein. Includes analysis of every individual
        # residue even if there is not an interaction occuring.
        # atom_num - include the number of the interacting atom - this option will significantly increase processing time
        # Atoms - which specific atoms of the ligand of interest 
        # reslabel - three letter code (ALA) or single letter code (A) for amino acids. 

        if atoms != None:
            atype = type(atoms) # Checks to see if atoms input is a list
            if atype != list:
                raise Exception(f'Atoms selection is not a valid input. Atoms should be a list of atoms, not {atype}')
        pdb = self.prot.pdb

        if smiles != None:
            lig_mol = Chem.rdmolfiles.MolFromPDBFile(self.ligand_mol)
            template = Chem.rdmolfiles.MolFromSmiles(smiles)
            lig_mol = AllChem.AssignBondOrdersFromTemplate(template, lig_mol)
        else:
            lig_mol = None

        chain_dfs = Protein_Analysis.check_multichain(pdb)
        interactions = []
        chain_num = 1
        for df in chain_dfs:
            df = df.assign(CHAIN = [chr(ord('@')+chain_num) for i in range(0, df.shape[0])])
            p_type = Protein_Analysis.identify_structure_type(df)
            if p_type == 'RNA' or p_type == 'DNA':
                proximity == True
            if reslabel == 'multiple':
                if len(chain_dfs) > 1:
                    df = df.assign(RESIDUE = df[['RES_NAME', 'RES_NUM', 'CHAIN']].agg('-'.join, axis=1))
                    chain_num += 1
                else:
                    df = df.assign(RESIDUE = df[['RES_NAME', 'RES_NUM']].agg('-'.join, axis=1))
            elif reslabel == 'single':
                res_labels = []
                if  p_type != 'RNA' and p_type != 'DNA':
                    amino_acids = {'ALA':'A','GLY':'G','VAL':'V','CYS':'C','PRO':'P','LEU':'L','ILE':'I','MET':'M','TRP':'W','PHE':'F',
                        'LYS':'K','ARG':'R','HIS':'H','SER':'S','THR':'T','TYR':'Y','ASN':'N','GLN':'Q','ASP':'D','GLU':'E'}
                    for i in df['RES_NAME'].values.tolist():
                        if i not in amino_acids:
                            res_labels.append(i)
                        else:
                            res_labels.append(amino_acids[i])
                else:
                    res_labels = df['RES_NAME']
                df['RES_SINGLE_NAME'] = res_labels
                if len(chain_dfs) > 1:
                    df = df.assign(RESIDUE = df[['RES_SINGLE_NAME', 'RES_NUM', 'CHAIN']].agg('-'.join, axis=1))
                else:
                    df = df.assign(RESIDUE = df[['RES_SINGLE_NAME', 'RES_NUM']].agg('-'.join, axis=1))
            else:
                raise Exception(f'{reslabel} is either not a valid weight option for residue labeling.')
            ligand = self.ligand.pdb                
    
            atom_cols = ['RES_NUM', 'ATOM_NUM', 'ATOM_TYPE', 'X', 'Y', 'Z']
            ratoms = df[atom_cols]
            residues = ratoms['RES_NUM'].unique()

            ligand = ligand.assign(ATOM_NAME_NUM = ligand[['ATOM_TYPE', 'ATOM_NUM']].agg('-'.join, axis=1))

            ligand_atoms = Protein_Analysis.find_atom(ligand, ligand['ATOM_NAME_NUM'].to_numpy(), atom_sel=atoms, ligand=True)

            max_dist = Protein_Analysis.max_diff(ligand_atoms)
            prot_atom_num = df['ATOM_NUM']
            lig_atom_num = ligand['ATOM_NUM']

            # End of data processing, start of classification of interactions

            for i in range(len(residues)):
                residue = df[df['RES_NUM']==residues[i]]
                residue_atoms = Protein_Analysis.find_atom(residue, residue['ATOM_TYPE'].to_numpy())
                for ligand_atom, residue_atom in itertools.product(ligand_atoms, residue_atoms):
                    dist = Protein_Analysis.euclidean_distance(ligand_atoms[ligand_atom], residue_atoms[residue_atom])
                    if dist > max_dist:
                        if atom_num == True:
                            lig_atom_num = ligand[(ligand['ATOM_TYPE'] == ligand_atom) & (ligand['X'] == ligand_atoms[ligand_atom][0])]['ATOM_NUM'].to_numpy()[0]
                            prot_atom_num = df[(df['ATOM_TYPE'] == residue_atom) & (df['X'] == residue_atoms[residue_atom][0])]['ATOM_NUM'].to_numpy()[0]
                            interactions.append([ligand_atom, lig_atom_num, residue_atom,residues[i],df['CHAIN_ID'].unique()[0],prot_atom_num, 'None', dist])
                        else:
                            interactions.append([ligand_atom, residue_atom,residue['RESIDUE'].unique()[0],df['CHAIN_ID'].unique()[0], 'None', dist])
                        break
                    interaction_type = Protein_Analysis.classify_interactions(dist, residue_atom, ligand_atom, residue, ligand, proximity, cutoff=cutoff, smiles=smiles, mol=lig_mol)
                    if atom_num == True:
                        lig_atom_num = ligand[(ligand['ATOM_TYPE'] == ligand_atom) & (ligand['X'] == ligand_atoms[ligand_atom][0])]['ATOM_NUM'].to_numpy()[0]
                        prot_atom_num = df[(df['ATOM_TYPE'] == residue_atom) & (df['X'] == residue_atoms[residue_atom][0])]['ATOM_NUM'].to_numpy()[0]
                        interactions.append([ligand_atom, lig_atom_num, residue_atom,residues[i],df['CHAIN_ID'].unique()[0],prot_atom_num, interaction_type, dist])
                    else: 
                        interactions.append([ligand_atom, residue_atom,residue['RESIDUE'].unique()[0], df['CHAIN_ID'].unique()[0], interaction_type, dist])

        df = pd.DataFrame(interactions)
        df = df.where(pd.notnull(df), 'None')
        
        if len(df.columns) == 7:
            df.columns = ['LIGAND_ATOM', 'LIG_ATOM_NUM', 'PROTEIN_ATOM', 'RESIDUE', 'CHAIN_ID','PROTEIN_ATOM_NUM', 'INTERACTION', 'DISTANCE']
        else: 
            df.columns = ['LIGAND_ATOM', 'PROTEIN_ATOM', 'RESIDUE', 'CHAIN_ID', 'INTERACTION', 'DISTANCE']

        if include_nones == True:
            df = df.reset_index(drop=True)
            return df
        else:
            df = df[df['INTERACTION'] != 'None']
            df = df.reset_index(drop=True)
            return df
        
    def check_multichain(prot):
        chains = prot['CHAIN_ID'].unique()
        dfs = []
        for i in chains:
            temp_df = prot[prot['CHAIN_ID'] == i]
            dfs.append(temp_df)
        return dfs
            
    def identify_structure_type(prot):
        amino_acids = {'ALA':'A','GLY':'G','VAL':'V','CYS':'C','PRO':'P','LEU':'L','ILE':'I','MET':'M','TRP':'W','PHE':'F',
                'LYS':'K','ARG':'R','HIS':'H','SER':'S','THR':'T','TYR':'Y','ASN':'N','GLN':'Q','ASP':'D','GLU':'E'}
        for i in prot['RES_NAME'].values.tolist():
            if i == 'C' or i == 'G':
                continue
            if i in amino_acids:
                return 'protein'
            else:
                if 'U' in prot['RES_NAME'].values.tolist():
                    return 'RNA'
                if 'T' in prot['RES_NAME'].values.tolist():
                    return 'DNA'
    
    def max_diff(atoms_coors):
        coors = list(atoms_coors.values())
        max_ = pdist(coors)
        return max_.max() + 7.5 + 5
        
    
    def clean_fingerprint(self, include_nones=False, proximity=False, weigh_by='interaction', atoms=None, reslabel='multiple', cutoff=5, smiles=None):
        df = Protein_Analysis.fingerprinting(self, include_nones=include_nones, proximity=proximity, atoms=atoms, reslabel=reslabel, cutoff=cutoff, smiles=smiles)
        interacting_residues = df['RESIDUE'].unique()
        fingerprint = []
        for i in interacting_residues:
            res_df = df[df['RESIDUE'] == i]
            distances = res_df['DISTANCE'].to_numpy()
            dist_sorted = np.sort(distances)
            res = res_df.loc[res_df['DISTANCE'] == dist_sorted[0]]
            if weigh_by == 'distance':
                fingerprint.append([res.iloc[0]['RESIDUE'], res.iloc[0]['INTERACTION'], res.iloc[0]['DISTANCE']])
            elif weigh_by == 'interaction':
                interactions = res_df['INTERACTION'].to_numpy()

                if 'Polar Charged' in interactions:
                    int_type = 'Polar Charged'
                elif 'H-Bond' in interactions:
                    int_type = 'H-Bond'
                elif 'Polar Uncharged' in interactions:
                    int_type = 'Polar Uncharged'
                elif 'Pi-stacking' in interactions:
                    int_type = 'Pi-stacking'
                elif 'Hydrophobic' in interactions:
                    int_type = 'Hydrophobic'
                elif 'Backbone' in interactions:
                    int_type = 'Backbone'
                elif 'None' in interactions:
                    int_type = 'None'
                elif 'Proximity' in interactions:
                    int_type = 'Proximity'
                
                res_int = res_df.loc[res_df['INTERACTION'] == int_type]
                fingerprint.append([res_int.iloc[0]['RESIDUE'], res_int.iloc[0]['INTERACTION'], res_int.iloc[0]['DISTANCE']])
            else:
                raise Exception(f'{weigh_by} is not a valid weight option for fingerprinting output. Please use interaction or distance')
                
        df = pd.DataFrame(fingerprint)
        
        df.columns = ['RESIDUE', 'INTERACTION', 'DISTANCE']

        if include_nones == False:
            df = df[df['INTERACTION'] != 'None']
        
        df = df.reset_index(drop=True)
                
        return df
        
    
    def euclidean_distance(p1, p2):
        dist = (((p1[0]-p2[0])**2)+((p1[1]-p2[1])**2)+((p1[2]-p2[2])**2))**(1/2)
        return dist
    
    def classify_interactions(distance, residue_atom, ligand_atom, residue, ligand, proximity, cutoff=5, smiles=None,  mol=None):
        lig_atom_list = ligand['ATOM_NAME_NUM'].to_list()
        lig_atom_list = [ x for x in lig_atom_list if "H" not in x ]
        if 'H' not in ligand_atom:
            ind = lig_atom_list.index(ligand_atom)

        ligand_atom = ligand_atom.split('-')[0]

        # print(lig_atom_list)
        # print(ind)
        # print(mol.GetAtomWithIdx(ind).GetSymbol())
        # print(Chem.MolToMolBlock(m2))

        if distance < cutoff:
            backbone_atoms = ['CA','C', 'N', 'O', 'HA']
            hydrophobic_residues = ['ALA', 'GLY', 'VAL', 'ILE', 'LEU', 'MET']
            pos_charged_residues = ['ARG', 'LYS', 'HIS']
            neg_charged_residues = ['GLU', 'ASP']
            aromatic_residues = ['TRP', 'TYR', 'PHE']
            # hbond_donors = ['HH', 'HE', 'HZ', 'HG', 'HO', 'HN']
            # hbond_acceptors = ['OD', 'OG']
            # aromatic_atoms = ['CE', 'CG', 'CD', 'CZ']
            charged_atoms = ['N', 'O', 'F', 'CL']
            # print(residue.iloc[0]['RES_NAME'])
            if smiles != None and 'H' not in ligand_atom:
                if (residue.iloc[0]['RES_NAME'] in aromatic_residues and (mol.GetAtomWithIdx(ind).GetIsAromatic() == True or ligand_atom in charged_atoms)) or (mol.GetAtomWithIdx(ind).GetIsAromatic() == True and (residue_atom in charged_atoms and (residue.iloc[0]['RES_NAME'] in pos_charged_residues or residue.iloc[0]['RES_NAME'] in neg_charged_residues))):
                    return 'Pi-stacking'
            if residue_atom in backbone_atoms:
                return 'Backbone'
            elif (('O' in residue_atom or 'N' in residue_atom or 'OP' in residue_atom) and 'C' not in ligand_atom) or (('O' in ligand_atom or 'N' in ligand_atom or 'OP' in ligand_atom) and 'C' not in residue_atom):
                if (residue.iloc[0]['RES_NAME'] in neg_charged_residues and 'N' in ligand_atom) or (residue.iloc[0]['RES_NAME'] in pos_charged_residues and 'O' in ligand_atom):
                    return 'Polar Charged'
                elif ('H' in residue_atom and ('O' in ligand_atom or 'N' in ligand_atom or 'OP' in ligand_atom)) or (('O' in residue_atom or 'N' in residue_atom or 'OP' in residue_atom) and ('H' in ligand_atom)):
                    return 'H-Bond'
                else:
                    return 'Polar Uncharged'
            elif ('C' in residue_atom or residue.iloc[0]['RES_NAME'] in hydrophobic_residues) and 'C' in ligand_atom:
                return 'Hydrophobic'
            elif proximity == True:
                return 'Proximity'
            else:
                return 'None'
        else:
            return 'None'

class Initialize_PDB:
    
    def __init__(self,pep):
        self.pep = pep
        self.get_from_pdb()

    def protein_info(self):
        table = self.pdb

        table['RES_NUM'] = table['RES_NUM'].astype(int)
        table['RES_NUM'] = table['RES_NUM'].astype(str)
        table['Residue_Num'] = table[['RES_NAME', 'RES_NUM']].agg('-'.join, axis=1)
        
        res = table['RES_NUM'].unique()
        
        code = table.loc[table['ATOM_TYPE'] == 'CA']
        three_code = code['RES_NAME'].to_numpy()
        
        atom_cols = ['Residue_Num', 'ATOM_NUM', 'ATOM_TYPE', 'X', 'Y', 'Z']
        
        atoms = table[atom_cols]

        length = len(res)
        amino_acids = {'ALA':'A','GLY':'G','VAL':'V','CYS':'C','PRO':'P','LEU':'L','ILE':'I','MET':'M','TRP':'W','PHE':'F',
                        'LYS':'K','ARG':'R','HIS':'H','SER':'S','THR':'T','TYR':'Y','ASN':'N','GLN':'Q','ASP':'D','GLU':'E'}
        seq_list = []
        for a in three_code:
            seq_list.append(amino_acids[a])
        seq = ''.join(seq_list)
        
        self.res_num = res
        self.three_code = three_code
        self.length = length
        self.seq = seq
        self.atoms = atoms
    
    def get_from_pdb(self):
        file = self.pep
        atom = []
        atom_num = []
        atom_type = []
        res_name = []
        chain_ID = []
        res_num = []
        x = []
        y = []
        z = []
        occupancy = []
        bfactor = []
        atom_name = []
        with open(file, "rt", newline='') as input_file:
            for line in input_file:
                if 'ATOM' in line or 'HETATM' in line:
                    atom.append(line[0:5].strip())
                    atom_num.append(line[6:11].strip())
                    atom_type.append(line[12:16].strip())
                    res_name.append(line[16:20].strip())
                    chain_ID.append(line[20:22].strip())
                    res_num.append(line[22:26].strip())
                    x.append(float(line[27:38].strip()))
                    y.append(float(line[38:46].strip()))
                    z.append(float(line[46:56].strip()))
                    occupancy.append(line[54:60].strip())
                    bfactor.append(line[60:68].strip())
                    atom_name.append(line[75:80].strip())
        df = pd.DataFrame([atom,atom_num,atom_type,res_name,chain_ID,res_num,x,y,z,occupancy,bfactor,atom_name]).transpose()
        df.columns = ['ATOM', 'ATOM_NUM', 'ATOM_TYPE', 'RES_NAME', 'CHAIN_ID', 'RES_NUM', 'X', 'Y', 'Z', 'OCCUPANCY', 'BFACTOR', 'ATOM_NAME']
        self.pdb = df


class Multipose_Docking_Fingerprinting:

    def __init__(self, prot, ligands):
        self.prot = prot
        self.ligands = ligands
        
    def docked_fingerprinting(self, proximity=False, weigh_by='interaction', atoms=None, reslabel='multiple', include_nones=False, cutoff=5, smiles=None):

        fingerprints = []

        for l in self.ligands:
            analysis = Protein_Analysis(self.prot, l)
            if weigh_by == 'interaction':
                fingerprint_data = analysis.clean_fingerprint(proximity=False, include_nones=True, atoms=atoms, reslabel=reslabel, cutoff=cutoff, smiles=smiles)
            elif weigh_by == 'distance':
                fingerprint_data = analysis.clean_fingerprint(proximity=False, include_nones=True, atoms=atoms, reslabel=reslabel, cutoff=cutoff, smiles=smiles)
            fingerprints.append(fingerprint_data)
        
        if len(fingerprints) == 0:
            print('File Empty - check to make sure output was correctly written')
            df = pd.DataFrame()
            return df
        df = fingerprints[0]

        df.columns = ['RESIDUE', 'INTERACTION POSE1', 'DISTANCE POSE1']
        range_tot = len(self.ligands)
        for i in range(1,range_tot):
            pose_num = i+1
            df[f'INTERACTION POSE{pose_num}'] = fingerprints[i]['INTERACTION']
            df[f'DISTANCE POSE{pose_num}'] = fingerprints[i]['DISTANCE']
        pose_cols = [f'INTERACTION POSE{i}' for i in range(1,range_tot+1)]
        
        df_update = df.copy()
        df_removenone = df_update[pose_cols].transpose()

        for i in df_removenone.columns:
            unique = df_removenone[i].unique()
            if len(unique) == 1 and unique[0] == 'None':
                df_update.drop([i], axis=0, inplace=True)
        
        if include_nones == False:
            df_update = df_update.reset_index(drop=True)
            return df_update
        elif include_nones == True:
            df = df.reset_index(drop=True)
            return df
        
    def calc_frequency(df):
        frequency = {}
        ints = [i for i in df.columns if 'INTERACTION' in i]

        cols = df['RESIDUE'].to_numpy()
        
        ints_df = df[ints].transpose()
        ints_df.columns = cols
        for i in cols:
            interactions = 0
            for v in ints_df[i].to_numpy():
                if v != 'None' and v != 'NaN':
                    interactions += 1
            frequency[i] = interactions/(len(ints)+1)

        return pd.DataFrame(frequency, index=[0]) 
    
    def heatmap(df, lig_labels, ax, color='rainbow'):
        col_names = df.columns
        _max = None
        for i in range(0,2):
            for num in df.max(axis=i):
                if _max is None:
                    _max = num
                else:
                    if num > _max:
                        _max = num

        graph_data = df.to_numpy()
        premask = np.tile(np.arange(graph_data.shape[1]), graph_data.shape[0]).reshape(graph_data.shape)

        images = []
        if color == 'rainbow' or color == 'interaction' or color == 'charged':
            for i in range(graph_data.shape[1]):
                col = np.ma.array(graph_data, mask = premask != i)
                if color == 'rainbow':    
                    im = ax.imshow(col, vmin=0, vmax=_max, cmap=Multipose_Docking_Fingerprinting.get_rainbow_colormap(i, len(col_names)))
                else:
                    im = ax.imshow(col, vmin=0, vmax=_max, cmap=Multipose_Docking_Fingerprinting.get_interactions_colormap(i, col_names[i], int_type=color))
                images.append(im)
        else:
            ax.imshow(graph_data, vmin=0, vmax=_max, cmap=color)
        ax.set(xticks=np.arange(len(df.columns)), yticks=np.arange(len(df)), xticklabels=col_names, yticklabels=lig_labels)
        plt.tight_layout()
        return ax
        
    def cross_analyze(*args):
        if type(args[0]) == list:
            args = args[0]
        df_freq = Multipose_Docking_Fingerprinting.calc_frequency(args[0])
        for arg in range(1,len(args)):
            df = Multipose_Docking_Fingerprinting.calc_frequency(args[arg])
            df_freq = pd.concat([df_freq, df])
        for i in df_freq.columns:
            unique = df[i].unique()
            if len(unique) == 1 and unique[0] == 0:
                df_freq.drop([i], axis=1, inplace=True)
        return df_freq
    
    def cluster(self):

        lig_COM = []

        ligs = self.ligands
        for l in ligs:
            analysis = Protein_Analysis(self.prot, l)
            lig = analysis.ligand.pdb
            lig = lig.assign(ATOM_NAME_NUM = lig[['ATOM_TYPE', 'ATOM_NUM']].agg('-'.join, axis=1))
            lig_atoms = Protein_Analysis.find_atom(lig, lig['ATOM_NAME_NUM'].to_numpy(), ligand=True)
            x_coors = []
            y_coors = []
            z_coors = []
            for key, value in lig_atoms.items():
                x_coors.append(value[0])
                y_coors.append(value[1])
                z_coors.append(value[2])
            lig_COM.append([sum(x_coors)/len(x_coors), sum(y_coors)/len(y_coors), sum(z_coors)/len(z_coors)])

        data = np.array(lig_COM)

        sil_score_max = -1 #this is the minimum possible score

        for n_clusters in range(2,9):
            model = KMeans(n_clusters = n_clusters, init='k-means++', max_iter=100, n_init=1)
            labels = model.fit_predict(data)
            sil_score = silhouette_score(data, labels)
            if sil_score > sil_score_max:
                sil_score_max = sil_score
                best_n_clusters = n_clusters

        model = KMeans(n_clusters=best_n_clusters)
        model.fit_predict(data)
        pred = model.fit_predict(data)

        return Counter(pred).keys(), Counter(pred).values()


    # Things that will eventually get moved to something else, just here for convenience for now
    
    def bar_graph(df, ax=None, title=None):
        if ax == None:
            fig, ax = plt.subplots(figsize=(3,3))
        df = df.columns.to_frame().T.append(df, ignore_index=True).transpose()
        df.columns = ['RESIDUES', 'FREQUENCY']
        ax.bar(df['RESIDUES'], df['FREQUENCY'], edgecolor='black', color='white')
        ax.tick_params(labelrotation=90, labelsize=10, axis='x')
        ax.tick_params(labelsize=10, axis='y') 
        ax.set_ylabel("Interaction Frequency", fontsize=10, fontweight="bold")
        ax.set_xlabel("Residue", fontsize=10, fontweight="bold")
        ax.set_ylim(0,1)
        if title != None:
            ax.set_title(title, weight='bold')
        return ax
    
    def get_rainbow_colormap(i, col_vals):
        col_range = [(255,255,255)]
        color_step = (i+1)/col_vals
        cmap = matplotlib.cm.get_cmap('jet')
        cmap_col = list(cmap(color_step))
        cmap_col.pop(3)
        set_col = []
        for i in cmap_col:
            col = int(round(i * 255))
            set_col.append(col)
        col_range.append(tuple(set_col))
        my_cmap = Multipose_Docking_Fingerprinting.create_colormap(col_range, bit=True)
        return my_cmap
    
    def get_interactions_colormap(i, col_names, int_type='interaction'):
        interaction_cols = {'charged': (125,114,178), 'polar': (222,131,0), 'aromatic': (240,236,66), 'hydrophobic': (0,117,115), 'positive' : (0,114,178), 'negative': (239,21,21), 'neutral': (103,113,115)} #blue 0,114,178
        interaction_type = {'charged': ['ARG', 'LYS', 'HIS', 'GLU', 'ASP', 'MG'], 'polar': ['SER', 'THR', 'ASN', 'GLN'], 'aromatic': ['TRP', 'TYR', 'PHE'], 'hydrophobic': ['ALA', 'GLY', 'VAL', 'ILE', 'LEU', 'MET', 'PRO', 'CYS']}
        amino_acids = {'ALA':'A','GLY':'G','VAL':'V','CYS':'C','PRO':'P','LEU':'L','ILE':'I','MET':'M','TRP':'W','PHE':'F',
                        'LYS':'K','ARG':'R','HIS':'H','SER':'S','THR':'T','TYR':'Y','ASN':'N','GLN':'Q','ASP':'D','GLU':'E', 'MG':'MG'}
        pos = ['ARG', 'LYS', 'HIS', 'MG']
        neg = ['GLU', 'ASP']
        aa = col_names.split('-')
        for key, int in interaction_type.items():
            for i in int:
                if len(aa[0]) == 1:
                    aa_int = amino_acids[i]
                else:
                    aa_int = i
                if aa_int == aa[0]:
                    if int_type == 'interaction':
                        interaction = key
                    else:
                        if i in pos:
                            interaction = 'positive'
                        elif i in neg:
                            interaction = 'negative'
                        else:
                            interaction = 'neutral'
        col_range = [(255,255,255)]
        col_range.append(interaction_cols[interaction])
        my_cmap = Multipose_Docking_Fingerprinting.create_colormap(col_range, bit=True)
        return my_cmap
    
    def create_colormap(colors, position=None, bit=False, reverse=False, name='custom_colormap'):
        from matplotlib.colors import LinearSegmentedColormap
        if not isinstance(colors, np.ndarray):
            colors = np.array(colors, dtype='f')
        if reverse:
            colors = colors[::-1]
        if position is not None and not isinstance(position, np.ndarray):
            position = np.array(position)
        elif position is None:
            position = np.linspace(0, 1, colors.shape[0])
        else:
            if position.size != colors.shape[0]:
                raise ValueError("position length must be the same as colors")
            elif not np.isclose(position[0], 0) and not np.isclose(position[-1], 1):
                raise ValueError("position must start with 0 and end with 1")
        if bit:
            colors[:] = [tuple(map(lambda x: x / 255., color)) for color in colors]
        cdict = {'red':[], 'green':[], 'blue':[]}
        for pos, color in zip(position, colors):
            cdict['red'].append((pos, color[0], color[0]))
            cdict['green'].append((pos, color[1], color[1]))
            cdict['blue'].append((pos, color[2], color[2]))
        return LinearSegmentedColormap(name, cdict, 256)
    
