# Written by Amanda K. Sharp
# In Southern Research
#
# Version 1
# Last updated 3/19/2024


import os
import numpy as np

def main():

    ligand = 'Target8.pdb' # Change this to your chemical target

    # The file below can be generated from my RCSB webscrapping script. If you
    # don't use that script, its just like a normal CSV file with the first column
    # representing the PDB IDs to use, and second column the gene name of the protein
    
    with open('scrubbed_pdbs.txt') as file: 
       
       for line in file:  
          
          pdb = line.split(',')[0]
          gene = line.split(',')[1]

          center, size = find_centroid_coors(pdb)

          gnina(f'{pdb}.pdb', ligand, center, size)

def gnina(receptor, ligand, center, size):
      
  # This function runs GNINA
  
  center_x=center[0]
  center_y=center[1]
  center_z=center[2]
  size_x=size[0]
  size_y=size[1]
  size_z=size[2]

  recept_name = receptor.split('.')[0]
  lig_name = ligand.split('.')[0]

  path_to_output = f'{recept_name}_{lig_name}_output.pdb'

  if os.path.isfile(path_to_output) == False:

    os.system(f'./../GNINA/gnina -r {receptor} -l {ligand} --center_x {center_x} --center_y {center_y} --center_z {center_z} --size_x {size_x} --size_y {size_y} --size_z {size_z} -o {path_to_output} --addH off --stripH off')

def find_centroid_coors(pdb):
  
  # This function probes all atom coordinates within a file and finds the outmost
  # atoms in the XYZ coordinates and the center to generate box size and center
  # to input into docking 
  
  coors = []
  with open(f'{pdb}.pdb', 'r') as file:
    for line in file:
      if line.startswith('HETATM') or line.startswith('ATOM'):
        l = line.split()
        coor = [float(i) for i in l[6:9]]
        coors.append(coor)
  coors_array = np.array([np.array(coor) for coor in coors])
  centroid = np.mean(coors_array[:,-3:], axis=0)

  centroid = [c for c in centroid]

  x_ = max([abs(x[0]-centroid[0]) for x in coors_array])*2
  y_ = max([abs(y[1]-centroid[1]) for y in coors_array])*2
  z_ = max([abs(z[2]-centroid[2]) for z in coors_array])*2

  size = [x_, y_, z_]

  return centroid, size

    
if __name__ == "__main__":
    main()