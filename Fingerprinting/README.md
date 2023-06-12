Fingerprinting.py
This script was written to generate data to analyze protein-ligand interactions. 
There are a few library dependencies required for this program to run. If you have Anaconda, you likely already have all these. Dependencies required:<br>
•	Numpy<br>
•	Pandas<br>
•	Itertools<br>
•	Matplotlib<br>
•	Scipy<br>

This script offers three different class methods for analyzing protein-ligand interactions: <br>
[fingerprinting.*Initialize_PDB*()](#fingerprinting.*Initialize_PDB*(input))<br>
[fingerprinting.*Protein_Analysis*()](#fingerprinting.*Protein_Analysis*(input))<br>
[fingerprinting.*Multipose_Docking_Fingerprinting*()](#fingerprinting.*Multipose_Docking_Fingerprinting*(input))<br>





# <code>fingerprinting.*Initialize_PDB*(input)</code>


##### Parameters: <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **input : *str*** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Input file name. PDB format required. Can be protein or ligand. <br>

##### Attributes: <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **pdb : *pandas dataframe*** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Holds raw information about input collected from PDB. <br>

#### Methods:

#### <code>protein_info()</code>
##### Attributes:<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **res_num : *list***<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Holds residue numbers for protein input.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **three_code : *list***<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Holds residue name with three letter code for protein input<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **length : *int***<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The number of residues in the protein input<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **seq : *str***<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Amino acid sequence for protein input<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **atoms : *pandas dataframe***<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Holds processed information from the protein input. <br>

##### Example:
```ruby
    protein = fingerprinting.Initialize_PDB('protein.pdb')
    protein.protein_info()
    print(protein.atoms)
```

##### Output: <br>
```ruby
     Residue_Num   ATOM_NUM   ATOM_TYPE      X       Y       Z
0          MET-1          1           N  22.072  10.521  25.004
1          MET-1          2          CA  22.547   9.419  24.101
2          MET-1          3           C  23.192  10.115  22.878
3          MET-1          4           O  22.763  11.195  22.464
4          MET-1          5          CB  21.363   8.482  23.747
...          ...        ...         ...     ...     ...     ...
3299     ASN-207       3300         1HB  38.515 -40.817 -16.251
3300     ASN-207       3301         2HB  38.882 -42.462 -16.703
3301     ASN-207       3302        1HD2  41.605   -40.7 -14.574
3302     ASN-207       3303        2HD2  39.874  -40.49 -14.521
3303      MG-203       3305          MG  15.565 -13.958   1.832
[3304 rows x 6 columns]
```


# <code>fingerprinting.*Protein_Analysis*(protein, ligand)</code>

This class provides methods to perform fingerprinting for a protein and a ligand. Protein and ligand must be in separate pdb files. 

##### Parameters: <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **protein : *str***<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Protein/receptor input file name. PDB format required.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **ligand : *str***<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Ligand input file name. PDB format required. Can be small molecule, peptide, or protein.<br>

#### Methods:

#### <code>fingerprinting(proximity=False, include_nones=False, atom_num=False, atoms=None, reslabel='multiple')</code>

Calculate the interactions between a protein (or biomolecule, DNA and RNA are supported) and a ligand (can be a small molecule, peptide, or protein). <br>

##### Parameters: <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **proximity : *bool, default=False*** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Include data that has no specific interaction classification but within 5 angstrom distance cutoff (may be interactions that would be considered repulsions. <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **include_nones : *bool, default=False*** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Returns data with all distances including "None" interactions. <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **atom_num : *bool, default=False*** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Returns data with atom numbers responsible for interactions. <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **atoms : *list, default=None*** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; A list of atoms of interest on the ligand. Output will only be residues that interact with those atoms.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **reslabel : *str, default='multiple'*** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Returns residue nomenclature are multiple, or three, code naming. Use *'single'* for naming with single letter naming. <br>

##### Returns: <br>

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **C : *pandas dataframe*** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Returns a pandas dataframe of all atomic interactions within a 5 angstrom cutoff. <br>
<br>

#### <code>clean_fingerprint(include_nones=False, proximity=False, weigh_by='interaction', atoms=None, reslabel='multiple')</code>
##### Parameters: <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **proximity : *bool, default=False*** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Include data that has no specific interaction classification but within 5 angstrom distance cutoff (may be interactions that would be considered repulsions. <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **include_nones : *bool, default=False*** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Returns data with all distances including "None" interactions. <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **weigh_by : *str, default='interaction'*** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; How the ranking of interaction is classified. Interaction selection will weigh interaction type for classification. Use *'distance'* for weighing classification of closest atom. <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **atoms : *list, default=None*** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; A list of atoms of interest on the ligand. Output will only be residues that interact with those atoms.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **reslabel : *str, default='multiple'*** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Returns residue nomenclature are multiple, or three, code naming. Use *'single'* for naming with single letter naming. <br>

##### Returns: <br>

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **C : *pandas dataframe*** <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Returns a pandas dataframe only of one interaction per residue within a 5 angstrom cutoff. <br>
<br>
