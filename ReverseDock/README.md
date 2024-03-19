This python script utilizes GNINA (required) to perform reverse docking. This will also identify the box size and center 
automatically, so no need to use external tools. At the moment, it covers the entirety of the protein, in future I will 
be optimizing and iteratively redocking for clustered hits with docking scores < -8 kcal/mol. 

To run this program, just use python reverse_dock.py. You will need to edit the contents of the file to your respective
chemical target and a text file (formatted like a csv file - comma deliminated) with the first column being the PDB ID
and the second column the gene name of the protein.
