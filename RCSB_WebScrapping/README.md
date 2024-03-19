This script is used to scrub a list of uniprot IDs to download PDB files. To run this script, use python get_pdbs.py. 
You will need to edit the main function of the script. Your input_file should be a text file with a list of line separated
uniprot IDs. You may also provide a ligand ID in case you only want to parse through a specific binding region of the list 
of proteins. Otherwise, remove the line that has ligand = "ATP" and the option screen_ligand=ligand from find_pdbs(). 

Unfortunately, Uniprot is not always super consistent with formatting. This script has been used on Kinases only at the moment. 
Updated will come in the future to hopefully streamline binding site identification. This script will also skip any protein
that does not have a PDB. I am working towards incorporating alphafold modeling for proteins that do not have a experimentally 
solved structure. 
