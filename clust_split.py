import numpy as np
import pandas as pd
import os

# Cluster PDB splitter

def split_clusters(cwd):
    files = []
    for (dirpath, dirnames, filenames) in os.walk(cwd):
        files.extend(filenames)
    for f in files:
        if 'pdb' in f: 
            with open(f, 'r') as infile:
                reader = infile.read()
                for i,part in enumerate(reader.split('MODEL')):
                    if i == 0:
                        print('Header for Top Cluster:\n' + part)
                    else:
                        with open('Clust' + str(i) + '_' + f, mode='w') as newfile:
                            newfile.write(part)
            infile.close()

if __name__ == '__main__':
    cwd = os.getcwd()
    split_clusters(cwd)