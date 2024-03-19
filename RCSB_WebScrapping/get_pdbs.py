# Written by Amanda K. Sharp
# In Southern Research
#
# Version 1
# Last updated 3/19/2024

import json
import requests
import urllib.parse
import requests
from pypdb import *
from pypdb.clients.pdb import pdb_client
import os

def main():

    uniprot_ids = []

    input_file = 'Uniprot_IDs.txt' 

    ligand = 'ATP'

    with open(input_file) as file:
        for line in file:
            up_ID = line.strip('\n')
            uniprot_ids.append(up_ID)

    pdbs, none, binding = find_pdbs(uniprot_ids, screen_ligand=ligand)

    if len(none) != 0:

        print(f'{none[0][0]} had no pdb files found that are dockable')

    get_pdbs(pdbs, binding)


def find_pdbs(uniprot_ids, screen_ligand=None):
    
    pdbs_reversedocking = []
    no_pdb = []
    binding = []

    for uniprot_id in uniprot_ids:

        print(f"Now collecting data for: {uniprot_id}")

        f = requests.get(f'https://rest.uniprot.org/uniprotkb/{uniprot_id}').json()

        gene = [f['genes'][0]['geneName']['value']]

        try:

            syn = [i['value'] for i in f['genes'][0]['synonyms']]

            genes = gene + syn

        except:

            genes = gene

        gene_name = f['genes'][0]['geneName']['value']

        print(f"{uniprot_id} identified as {gene_name}")

        ligands = []
        binding_site = []

        domains = [i['type'] for i in f['features']]

        if 'Binding site' in domains:

            bs_ind = domains.index('Binding site')

            if screen_ligand != None:

                if f['features'][bs_ind]['ligand']['name'] == screen_ligand:

                    binding_site.append(f['features'][bs_ind]['location']['start']['value'])
                    binding_site.append(f['features'][bs_ind]['location']['end']['value'])
                    ligands.append(f['features'][bs_ind]['ligand']['name'])
            else:
                binding_site.append(f['features'][bs_ind]['location']['start']['value'])
                binding_site.append(f['features'][bs_ind]['location']['end']['value'])
                ligands.append(f['features'][bs_ind]['ligand']['name'])

        # else:

        #     print(f'Binding site for {gene_name} was unable to be identified. You may need to mannually select a structure file for this gene.')
        #     break
        
        if len(binding_site) != 0:

            binding.append(binding_site)

            jsonfile = requests.get(f'https://rest.uniprot.org/uniprotkb/search?query={uniprot_id}&fields=structure_3d').json()

            pdb_ids = []

            for key, value in jsonfile.items():
                for val in value:
                    for i in val['uniProtKBCrossReferences']:
                        num_hits = 0
                        if ',' not in i['properties'][2]['value']:
                            try:
                                length_chains = i['properties'][2]['value'].split('=')[1].split('-')
                                for residues in binding_site:
                                    if int(length_chains[0]) <= residues <= int(length_chains[1]):
                                        num_hits += 1
                                    else:
                                        continue
                                if len(binding_site) == num_hits:
                                    pdb_ids.append(i['id'])
                            except:
                                print(i['properties'][2]['value'])

            selected = select_pdb(pdb_ids)

            if selected != None:
                pdbs_reversedocking.append([selected[0], genes, selected[1]])

            else: 
                no_pdb.append([uniprot_id, genes])

    else:
        no_pdb.append([uniprot_id, genes])

    return pdbs_reversedocking, no_pdb, binding




                    
def select_pdb(pdb_ids):

    ''' This function prioritizes the selection of a pdb file with a cocrystallized ligand '''

    if len(pdb_ids) != 0:

        resolutions = []

        for pdb_id in pdb_ids:
            pdb_info = get_all_info(pdb_id)

            try:

                ligands = [i['comp_id'] for i in pdb_info['rcsb_binding_affinity']]
                ligands = list(set(ligands))

            except:

                ligands = []  

            finally:
                if len(ligands) != 0:
                    return pdb_id, ligands
                
                else:
                    try:
                        resolutions.append([pdb_id, pdb_info['rcsb_entry_info']['resolution_combined'][0]])
                    
                    except:
                        pass

                    if pdb_id == pdb_ids[-1]:
                        try:
                            min_res = min(resolutions)
                            return min_res[0], None
                        
                        except:
                            return None 

        
def get_pdbs(pdbs, binding):

    amino_acids = ['ALA','GLY','VAL','CYS','PRO','LEU','ILE','MET','TRP','PHE','LYS','ARG','HIS',
                   'SER','THR','TYR','ASN','GLN','ASP','GLU', 'CYX', 'HID', 'AMET']


    filter_id = 'GENE'

    filters = ['GENE', 'ORGANISM_SCIENTIFIC']

    for p in pdbs:

        pdb = p[0]
        gene_list = p[1]
        pdb_file = pdb_client.get_pdb_file(pdb)
        ligand = p[2]

        pdb_lines = pdb_file.split('\n')
        keep_chain = []
        mols = 0

        try:

            for line in pdb_lines:

                if 'SOURCE' in line:
                    if 'MOL_ID' in line:
                        mols +=1

                    if filter_id in line:

                        filtered = line.split(':')[1].split(';')[0].strip()
            
                        if filtered in gene_list:
                            keep_chain.append(mols)

                        else:
                            continue
                
        except:

            pass

        if len(keep_chain) == 0:

            keep_chain = [1]


        chain_num = 1
        chains = []

        for line in pdb_lines:
            if line.startswith('HETATM') or line.startswith('ATOM'):
                chain_id = line.split()[4]
                if chain_id not in chains:
                    chains.append(chain_id)

        check = ['2BUJ', '3PLS', '5GRN', '7MON']

        with open(f'{pdb}_apo.pdb', 'a') as f:
            print(f'Writing PDB file for {pdb}')
            for line in pdb_lines:
                for keep in keep_chain:
                    if line is not None:
                        if line.startswith('HETATM') or line.startswith('ATOM'):
                            if chains[keep-1] == line.split()[4] and 'HOH' not in line: # 
                                    amino_acid = line.split()[3]
                                    if len(amino_acid) > 3:
                                        # print(amino_acid)
                                        amino_acid = amino_acid[-3]
                                    if amino_acid in amino_acids:
                                        # if ligand is not None and amino_acid in ligand:
                                        #     f.write(f'{line}\n')
                                        # else: 
                                        f.write(f'{line}\n')

        with open('scrubbed_pdbs.txt', 'a') as f:
            f.write(f'{pdb},{gene_list[0]}\n')


if __name__ == "__main__":
    main()