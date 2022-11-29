# Written by Amanda K. Sharp
# In the Brown Lab
#
# Version 2
# Last updated 3/18/21

import pandas as pd 
import numpy as np
pd.set_option('display.max_rows', None)


def All_Contact(df, num_poses):
    allcontact = df[df['interaction type'].str.contains('contact')]
    if len(num_poses) != 1:
        poses = [f'Pose{i}' for i in range(1,len(num_poses)+1)]
        poses = allcontact[poses]
        for i in range(1,len(num_poses)+1):
            poses[f'Pose{i}'] = poses[f'Pose{i}'].astype(int)
        allcontact['total'] = poses.sum(axis=1)
        allcontact['normalized'] = allcontact['total'].div(len(num_poses))
    else:
        allcontact['total'] = allcontact['Pose1']
        allcontact['normalized'] = allcontact['total']
    return(allcontact)

def Charged(df, num_poses):
    polarcharged = df[df['interaction type'].str.contains('charged')]
    if len(num_poses) != 1:
        poses = [f'Pose{i}' for i in range(1,len(num_poses)+1)]
        poses = polarcharged[poses]
        for i in range(1,len(num_poses)+1):
            poses[f'Pose{i}'] = poses[f'Pose{i}'].astype(int)
        polarcharged['total'] = poses.sum(axis=1)
        polarcharged['normalized'] = polarcharged['total'].div(len(num_poses)) 
    else:
        polarcharged['total'] = polarcharged['Pose1']
        polarcharged['normalized'] = polarcharged['total']
    return(polarcharged)
    
def Polar(df, num_poses):
    polar = df[df['interaction type'].str.contains('polar')]
    if len(num_poses) != 1:
        poses = [f'Pose{i}' for i in range(1,len(num_poses)+1)]
        poses = polar[poses]
        for i in range(1,len(num_poses)+1):
            poses[f'Pose{i}'] = poses[f'Pose{i}'].astype(int)
        polar['total'] = poses.sum(axis=1)
        polar['normalized'] = polar['total'].div(len(num_poses)) 
    else:
        polar['total'] = polar['Pose1']
        polar['normalized'] = polar['total']
    return(polar)
    
def Hydrophobic(df, num_poses):
    hydrophobic = df[df['interaction type'].str.contains('hydrophobic')]
    if len(num_poses) != 1:
        poses = [f'Pose{i}' for i in range(1,len(num_poses)+1)]
        poses = hydrophobic[poses]
        for i in range(1,len(num_poses)+1):
            poses[f'Pose{i}'] = poses[f'Pose{i}'].astype(int)
        hydrophobic['total'] = poses.sum(axis=1)
        hydrophobic['normalized'] = (hydrophobic['total']/len(num_poses))
    else:
        hydrophobic['total'] = hydrophobic['Pose1']
        hydrophobic['normalized'] = hydrophobic['total']
    return(hydrophobic)
    
def Backbone(df, num_poses):
    Backbone = df[df['interaction type'].str.contains('backbone')]
    if len(num_poses) != 1:
        poses = [f'Pose{i}' for i in range(1,len(num_poses)+1)]
        poses = Backbone[poses]
        for i in range(1,len(num_poses)+1):
            poses[f'Pose{i}'] = poses[f'Pose{i}'].astype(int)
        Backbone['total'] = poses.sum(axis=1)
        Backbone['normalized'] = (Backbone['total']/len(num_poses))
    else:
        Backbone['total'] = Backbone['Pose1']
        Backbone['normalized'] = Backbone['total']
    return(Backbone)
