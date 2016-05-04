"""
Reads in the measurements, adds a column for binary expression, pickles it,
and writes it to a csv.
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

    # get the codes from the names
    dict_df = pd.read_excel('dict.xlsx')
    code_dict = {n.lower():c for n,c in zip(dict_df['name'], dict_df['code'])}
    codes = [code_dict[n] for n in df['name']]
    df['code'] = codes

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

    # make columns for log_mKate and log_gfp
    df['log_mKate'] = np.log(df['mKate_mean'])
    df['log_GFP'] = np.log(df['GFP_mean'])

#     # calculate percentage of residues in common with each parent
#     parent_names = ['c1c2', 'cschrimson', 'cheriff']
#     for pn in parent_names:
#         parent_index = df['name'] == pn
#         parent = df[parent_index]['sequence']
#         df[pn+'_fraction'] = [identity(parent.item(),s) for s in df['sequence']]

    # make split for mkate_mean
    cut = df['mKate_mean'].fillna(0).median()
    bin_mKate = [1 if m > cut else -1 for m in df['mKate_mean']]
    df['bin_mKate_2'] = bin_mKate
    cut = df['mKate_mean'].dropna().median()
    bin_mKate = [np.nan if np.isnan(r) else\
                     1 if r > cut else -1 for r in df['mKate_mean']]
    df['bin_mKate'] = bin_mKate
    # make split for GFP
    cut = df['GFP_mean'].dropna().median()
    bin_GFP = [np.nan if np.isnan(r) else\
                     1 if r > cut else -1 for r in df['GFP_mean']]
    df['bin_GFP'] = bin_GFP
    # make split for sum_ratio
    cut = df['sum_ratio'].dropna().median()
    bin_sum_ratio = [np.nan if np.isnan(r) else\
                     1 if r > cut else -1 for r in df['sum_ratio']]
    df['bin_sum_ratio'] = bin_sum_ratio
    # make split for intensity_ratio
    cut = df['intensity_ratio'].dropna().median()
    bin_intensity_ratio = [np.nan if np.isnan(r) else\
                     1 if r > cut else -1 for r in df['intensity_ratio']]
    df['bin_intensity_ratio'] = bin_intensity_ratio

    cut = df[df['name']=='cheriff']['mKate_mean'].values
    above_parent = [1 if m >= cut else -1 for m in df['mKate_mean']]
    df['mKate_above_parent'] = above_parent

    cut = df[df['name']=='c1c2']['GFP_mean'].values
    above_parent = [1 if m >= cut else -1 for m in df['GFP_mean']]
    df['GFP_above_parent'] = above_parent

    cut = df[df['name']=='c1c2']['sum_ratio'].values
    above_parent = [1 if m >= cut else -1 for m in df['sum_ratio']]
    df['sum_ratio_above_parent'] = above_parent

    cut = df[df['name']=='c1c2']['intensity_ratio'].values
    above_parent = [1 if m >= cut else -1 for m in df['intensity_ratio']]
    df['intensity_ratio_above_parent'] = above_parent

    # pickle
    with open(os.path.join(dir, args.name.split('/')[0]
                           + '/props.pkl'), 'wb') as f:
        pickle.dump(df, f)
    # save as a csv
    with open(os.path.join(dir, args.name.split('/')[0]
                           + '/props.csv'), 'wb') as f:
        df.to_csv(f)

def identity(str1, str2):
    '''
    Find fraction of characters shared by the two strings
    '''
    if len(str1) != len(str2):
        raise ValueError ('Arguments must be strings of the same length')
    return sum([1.if a==b else 0. for a,b in zip(str1,str2)])/float(len(str1))

if __name__ == "__main__":
    main()
