### Fingerprint Analysis from Schrodinger-Maestro Fingerprinting output

These scripts were written to visualize the fingerprinting output from Schrodinger-Maestro (https://www.schrodinger.com/). Files should be saved using the nomenclature of <code>proteinname_ligandname_fingerprint.csv</code> for the python script to recognize. Each ligand should be saved into individual csv files. This can be all poses from docking in one csv files or top/best pose of the ligand in one csv file.

To identify which arguments are required to run this script, type: <code>python fingerprint.py -h</code> This will show you what arguments are required, optional, options that are accepted for each argument, and what is customizable in the script.

Bar graphs are helpful when looking at a specific ligand to see the interaction propensity for all the docked poses. This should only be used if using all docked poses. If you are wanting bar graphs for multiple ligands, we recommend only adding one csv in at first (only having once csv in the folder you are working in), figuring out the dimensions of the graph (see if the default fits well) using <code>-s yes</code>, then adding the rest in once dimensions have been worked out. An example for producing a bar graph using this methods is <code>python fingerprint.py -r protein.pdb -i all -g bar</code>.

Example Bar Graph:

<p align="center"><a href="https://imgur.com/snxi68o"><img src="https://i.imgur.com/snxi68o.png" title="source: imgur.com"  width="400" /></a></p>

A heatmap is useful if you have multiple ligands to show interactions for all the poses of the ligands (but top pose can also be used). If all poses from molecular docking is used, the darker/more vibrant color shows that more poses interact with that residue. The lighter the color, the fewer poses interact with that residue. To make a heatmap, use <code>python fingerprint.py -r protein.pdb -i all -g heatmap</code>. Again it is recommended to play around with the figure dimensions using the options in the script (-h will help you figure out exactly what those argument flags are). For this, <code>-x 10</code> and <code>-y 5 </code>is a reasonable size and adjust from there.

Example Heatmap:

<p align="center"><a href="https://imgur.com/db2cTYV"><img src="https://i.imgur.com/db2cTYV.png" title="source: imgur.com"  width="800" /></a></p>
