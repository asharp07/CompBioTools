# Written by Amanda K. Sharp
# In the Brown Lab
#
# Version 2
# Last updated 3/18/21


import renumber_fromPDB as RPDB
from custom_colormaps import create_colormap
import import_data
import argparse
import pandas as pd 
import numpy as np
import matplotlib 
import matplotlib.pyplot as plt
pd.options.mode.chained_assignment = None
import os


def main():
    residue_key = RPDB.Protein_Key(receptor, my_dir)
    res_match, ligands = import_data.gather_finger_data(residue_key, my_dir, interaction, datatype, ignorechain)
    data = np.concatenate(res_match)
    for i in range(len(res_match)-1):
        data = np.delete(data, i+2, 0)
    df = pd.DataFrame(data)
    transform_data(df, receptor, ligands)
    
def transform_data(data, protein, ligands):
    data.columns = data.iloc[0]
    data = data.reset_index()
    data = data.iloc[1: , :]
    data = data.astype(float)
    data = (data.loc[:, (data != 0).any(axis=0)])
    del data['index']
    col_names = data.columns
    graph_data = data.to_numpy()
    receptor, ex = protein.split('.pdb')
    try:
        graph_dict[graphtype](graph_data, col_names, receptor, ligands)
    except KeyError:
        print('Error in Graph Type')

def fingerprint_bar_graph(graph_data, col_names, receptor, ligands): 
    for i in range(len(graph_data)):
        cols = []
        ligand = ligands[i]
        for x in graph_data[i]:
            if x > 0.7:
                cols.append(color1)
            else:
                cols.append(color2)
        fig, ax = plt.subplots(figsize=(xlen, ylen))
        plt.bar(col_names, graph_data[i], color=cols, edgecolor='black')
        plt.xticks(rotation=90, fontsize=10)
        plt.yticks(fontsize=12) 
        plt.ylabel("Interaction Propensity", fontsize=14, fontweight="bold")
        plt.xlabel("Residue", fontsize=14, fontweight="bold")
        if threshold == 'yes':
            plt.axhline(y=0.7, color='gray', linestyle='--', zorder=0)
        plt.title(f'{receptor} and {ligand}', fontsize=18, fontweight="bold")
        plt.tight_layout()
        if show == 'yes':
            plt.show()
        else:
            plt.savefig(f'{receptor}_{ligand}_fingerprint_bargraph.png')

def fingerprint_heatmap(data, col_names, receptor, ligands):
    graph_data = pd.DataFrame(data)
    df = graph_data
    _max = None
    for i in range(0,2):
        for num in df.max(axis=i):
            if _max is None:
                _max = num
            else:
                if num > _max:
                    _max = num
    graph_data = graph_data.astype(float)
    graph_data = graph_data.to_numpy()
    premask = np.tile(np.arange(graph_data.shape[1]), graph_data.shape[0]).reshape(graph_data.shape)
    fig,ax = plt.subplots(figsize=(xlen, ylen))
    images = []
    for i in range(graph_data.shape[1]):
        col = np.ma.array(graph_data, mask = premask != i)
        im = ax.imshow(col, vmin=0, vmax=_max, cmap=get_colormap(i, len(col_names)))
        images.append(im)
    ax.set(xticks=np.arange(len(df.columns)), yticks=np.arange(len(df)), xticklabels=col_names, yticklabels=ligands)
    plt.xticks(rotation=90) 
    plt.ylabel("Ligand", fontsize=14, fontweight="bold")
    plt.xlabel("Residue", fontsize=14, fontweight="bold")
    plt.tight_layout()
    if show == 'yes':
        plt.show()
    else:
        plt.savefig(f'{receptor}_fingerprint_heatmap.png')

def get_colormap(i, col_vals):
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
    my_cmap = create_colormap(col_range, bit=True)
    return my_cmap

if __name__ == '__main__':
    interaction_type_options = ['all', 'charged', 'polar', 'hydrophobic', 'backbone']
    graph_type_options = ['heatmap', 'bar']
    parser = argparse.ArgumentParser(description = 'Fingerprinting Graph Customization Options')
    parser.add_argument('-i', '--interaction', help='(REQUIRED) Allowed interaction types are '+ ', '.join(interaction_type_options), required=True)
    parser.add_argument('-c1', '--color1', help='(OPTIONAL) High Interaction Color (hex) - Only required for bar graph', required=False, default='#6197AD')
    parser.add_argument('-c2', '--color2', help='(OPTIONAL) Low Interaction Color (hex) - Only required for bar graph', required=False, default='gray')
    parser.add_argument('-g', '--graphtype', help='(REQUIRED) Graph Type - Allowed graph types are '+', '.join(graph_type_options), required=True)
    parser.add_argument('-d', '--datatype', help='(OPTIONAL) Do you want your data counts to be total or normalized (0-1)? (total or normalized)', required=False, default='total')
    parser.add_argument('-r', '--receptor', help='(REQUIRED) Receptor file input', required=True)
    parser.add_argument('-x', '--xlen', help='(OPTIONAL) For figure dimensions - x length of the figure', required=False, default=6)
    parser.add_argument('-y', '--ylen', help='(OPTIONAL) For figure dimensions - y length of the figure', required=False, default=5)
    parser.add_argument('-t', '--threshold', help='(OPTIONAL) Add a threshold line to your graph (yes or no)', required=False, default='no')
    parser.add_argument('-s', '--show', help='(OPTIONAL) Display the graph in a popup rather than autosaving (yes or no)', required=False, default='no')
    parser.add_argument('-ch', '--ignorechain', help='(OPTIONAL) Do not include chain ID in the graph (not recommended for multichain protein complexes. Use yes if you want to ignore chain ID in graph)', required=False, default='no')

    args = parser.parse_args()
    interaction = args.interaction
    color1 = args.color1
    color2 = args.color2
    graphtype = args.graphtype
    datatype = args.datatype
    receptor = args.receptor
    xlen = args.xlen
    ylen = args.ylen
    threshold = args.threshold
    show = args.show
    ignorechain = args.ignorechain

    xlen = np.array(xlen, dtype=float)
    ylen = np.array(ylen, dtype=float)

    graph_dict = {'bar' : fingerprint_bar_graph, 'heatmap' : fingerprint_heatmap}

    my_dir = os.getcwd() 
    main()