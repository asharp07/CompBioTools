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


##### Example:
This analysis looks at residues interacting with the phosphate atoms of the ligand. <br> 
```ruby
    protein = 'protein.pdb'
    ligand = 'ligand.pdb'
    analysis = fingerprinting.Protein_Analysis(protein, ligand)  
    fingerprint = analysis.fingerprinting(atoms=['P'])
    
    print(fingerprint)
```

##### Output: <br>
```ruby
   LIGAND_ATOM PROTEIN_ATOM  RESIDUE CHAIN_ID      INTERACTION  DISTANCE
0          PA1           NE   ARG-20        A    Polar Charged  4.245834
1          PA1          NH1   ARG-20        A    Polar Charged  4.097387
2          PA1          NH2   ARG-20        A    Polar Charged  3.680664
3          PA2          NH1   ARG-20        A    Polar Charged  4.807463
4          PA2          NH2   ARG-20        A    Polar Charged  4.400783
...        ...          ...      ...       ...             ...       ...
41         PA5           HH  TYR-156        A         Aromatic  4.381776
42         PB5           CZ  TYR-156        A         Aromatic  4.808794
43         PB5           OH  TYR-156        A  Polar Uncharged  3.668233
44         PB5           HH  TYR-156        A         Aromatic  3.398084
[45 rows x 6 columns]
```

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


##### Example:
```ruby
    protein = 'protein.pdb'
    ligand = 'ligand.pdb'
    analysis = fingerprinting.Protein_Analysis(protein, ligand)  
    fingerprint = analysis.clean_fingerprint(atoms=['P'])
    
    print(fingerprint)
```

##### Output: <br>
```ruby
   RESIDUE      INTERACTION  DISTANCE
0   ARG-20    Polar Charged  4.245834
1   SER-54  Polar Uncharged  3.688160
2   ASN-56         Backbone  4.673043
3   VAL-61      Hydrophobic  4.111853
4   LYS-64    Polar Charged  3.328457
5   GLY-65      Hydrophobic  3.802741
6  PHE-101         Aromatic  4.930385
7  SER-103  Polar Uncharged  4.799216
8  ARG-140    Polar Charged  4.848395
9  TYR-156  Polar Uncharged  4.196831
[10 rows x 3 columns]
```

# <code>fingerprinting.*Multipose_Docking_Fingerprinting*(protein, ligands)</code>

This class provides methods to perform fingerprinting for a protein and a ligand. Protein and ligand must be in separate pdb files. 

##### Parameters: <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **protein : *str***<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Protein/receptor input file name. PDB format required.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **ligands : *list***<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; List of ligand input file names. PDB format required. Can be small molecule, peptide, or protein.<br>

#### Methods:

#### <code>docked_fingerprinting(proximity=False, weigh_by='interaction', atoms=None, reslabel='multiple', include_nones=False)</code>

Calculate the interactions between a protein (or biomolecule, DNA and RNA are supported) and a ligand (can be a small molecule, peptide, or protein). <br>

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
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Returns a pandas dataframe of all atomic interactions within a 5 angstrom cutoff. <br>
<br>

##### Example:
```ruby
    protein = 'protein.pdb'
    ligands = [f'ligand_pose{i}.pdb' for i in range(1, 10)]
    analysis = fingerprinting.Multipose_Docking_Fingerprinting(protein, ligands)
    df = analysis.docked_fingerprinting()
    print(df)
```

##### Output: <br>
```ruby
    RESIDUE INTERACTION POSE1  DISTANCE POSE1 INTERACTION POSE2  ...  INTERACTION POSE8 DISTANCE POSE8  INTERACTION POSE9 DISTANCE POSE9
0     ALA-6              None       12.211658       Hydrophobic  ...               None      15.763134        Hydrophobic       4.305532
1     ARG-7              None       10.445213          Backbone  ...               None      14.644728           Backbone       4.904095
2     THR-8              None        9.019689          Backbone  ...               None      14.533458           Backbone       3.993462
3     GLY-9              None        7.350499       Hydrophobic  ...               None      13.703370        Hydrophobic       4.970275
4    ARG-10              None        7.554643              None  ...               None      13.596932               None      14.608713
...     ...               ...             ...               ...  ...                ...            ...                ...            ...
26  GLU-136              None       16.099464              None  ...               None      18.222701      Polar Charged       4.872008
27  ARG-140     Polar Charged        4.848395     Polar Charged  ...      Polar Charged       4.550965      Polar Charged       4.736821
28  TYR-156   Polar Uncharged        4.196831   Polar Uncharged  ...    Polar Uncharged       4.838458               None      11.046769
29  TRP-158          Aromatic        4.997566   Polar Uncharged  ...    Polar Uncharged       4.780807               None      13.413603
[30 rows x 19 columns]
```
