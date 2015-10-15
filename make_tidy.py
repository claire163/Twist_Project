"""
Reads in the measurements, adds a column for binary expression, pickle it,
and write it to a csv.
"""

import pandas as pd
import argparse
import os
import chimera_tools
import pickle
from sklearn.cluster import KMeans
from sys import exit
import numpy as np


def main():
    dir = os.path.dirname(__file__)
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--name',required=True)

    args = parser.parse_args()

    df = pd.read_excel(args.name)

    # make names lowercase
    df['name'] = [s.lower() for s in df['name']]


    # zero index codes
    df['code'] = [chimera_tools.zero_index(c) for c in df['code']]

    # add column for binary expression
    # first, get indices where 'mkate_mean' is NaN
    expressed = pd.Series.notnull(df['mKate_mean'])
    # make the expression column
    df['expression'] = [1 if e else -1 for e in expressed]

    # generate sequence as string
    a_and_c = os.path.join(dir, 'alignment_and_contacts.pkl')
    sample_space, contacts = pickle.load(open(a_and_c))
    c_assignments_file = os.path.join(dir,'clibrary.output')
    n_assignments_file = os.path.join(dir,'nlibrary.output')
    c_assignments_dict = chimera_tools.load_assignments(c_assignments_file)
    n_assignments_dict = chimera_tools.load_assignments(n_assignments_file)
    dict_file = os.path.join(dir,'dict.xlsx')
    name_dict = chimera_tools.make_name_dict(dict_file)

    seqs = []
    for name in df['name']:
        if name[0] in set(['n', 'N']):
            seqs.append(''.join(chimera_tools.make_sequence(name_dict[name], n_assignments_dict, sample_space)))
        else:
            seqs.append(''.join(chimera_tools.make_sequence(name_dict[name], c_assignments_dict, sample_space)))
    df['sequence'] = seqs


    # calculate percentage of residues in common with each parent
    parent_names = ['c1c2', 'cschrimson', 'cheriff']
    for pn in parent_names:
        parent_index = df['name'] == pn
        parent = df[parent_index]['sequence']
        df[pn+'_fraction'] = [identity(parent.item(),s) for s in df['sequence']]


    # get the raw gfp values from the ratio and the mKate values
    df['gfp'] = df['sum_ratio']*df['mKate_mean']

    # cluster expressors
    kmeans = KMeans(n_clusters=3)
    exp = np.transpose(np.array(df[expressed]['mKate_mean']))
    exp = np.matrix(exp).T
    kmeans.fit(exp)
    centers = sorted(kmeans.cluster_centers_)
    print centers
    exp_level = [-1 if df.loc[i]['expression']==-1 else which_cluster(centers,df.loc[i]['mKate_mean'])\
                 for i in df.index]
    df['exp_level'] = exp_level

    # cluster for localization
    kmeans = KMeans(n_clusters=3)
    loc = np.array(df[expressed]['sum_ratio'])
    loc = np.matrix(loc).T
    kmeans.fit(loc)
    centers = sorted(kmeans.cluster_centers_)
    loc_level = [-1 if df.loc[i]['expression']==-1 else which_cluster(centers,df.loc[i]['sum_ratio'])\
                 for i in df.index]
    df['loc_level'] = loc_level
    print centers

    # pickle
    with open(os.path.join(dir, args.name.split('.xlsx')[0] + '.pkl'), 'wb') as f:
        pickle.dump(df, f)
    # save as a csv
    with open(os.path.join(dir, args.name.split('.xlsx')[0] + '.csv'), 'wb') as f:
        df.to_csv(f)


def which_cluster(centers, point):
    '''
    Return index of cluster closest to the point
    '''
    distances = [(c-point)**2 for c in centers]
    return np.argmin(distances)

def identity(str1, str2):
    '''
    Find fraction of characters shared by the two strings
    '''
    if len(str1) != len(str2):
        raise ValueError ('Arguments must be strings of the same length')
    return sum([1.if a==b else 0. for a,b in zip(str1,str2)])/float(len(str1))

if __name__ == "__main__":
    main()
