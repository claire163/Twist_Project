"""
Reads in the measurements, adds a column for binary expression, pickles it,
and writes it to a csv.
"""

import pandas as pd
import argparse
import os
import sys
sys.path.append('/Users/kevinyang/Documents/Projects/GPModel')
sys.path.append ('/Users/seinchin/Documents/Caltech/Arnold Lab/Programming tools/GPModel')
import chimera_tools
import pickle
from sklearn.cluster import KMeans
from sys import exit
import numpy as np

def split(df, column_name, how):
    try:
        if how == 'median':
            cut = df[column_name].dropna().median()
            bin_ = [np.nan if np.isnan(r) else
                    1 if r > cut else -1
                    for r in df[column_name]]
            df['bin_' + column_name.split('_mean')[0]] = bin_
        elif how == 'parents':
            parents = ['cschrimson', 'cheriff', 'c1c2']
            cut = min(df[df['name'].isin(parents)][column_name].values)
            bin_ = [1 if m >= cut else -1 for m in df[column_name]]
            df[column_name.split('_mean')[0] + '_above_parent'] = bin_
    except KeyError:
        pass

def make_subblocks(code, conversion_dict):
    subblocks = ''.join([code[conversion_dict[k]] for k in range(len(conversion_dict))])

def main():
    dir = os.path.dirname(__file__)
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--name',required=True)
    parser.add_argument('-c', '--min_cells', required=False, type=int)

    args = parser.parse_args()

    df = pd.read_excel(args.name)
    if args.min_cells is not None:
        df = df[df['cell'] > args.min_cells]
    # make names lowercase
    df['name'] = [s.lower() for s in df['name']]

    # get the codes from the names
    dict_df = pd.read_excel('dict.xlsx')
    code_dict = {n.lower():c for n,c in zip(dict_df['name'], dict_df['code'])}
    codes = [code_dict[n] for n in df['name']]
    df['code'] = codes

    # zero index codes
    df['code'] = [chimera_tools.zero_index(c) for c in df['code']]

    # generate subblocks
    sub_df = pd.read_csv('df_s_c_n_sblock.csv')
    for c in ['s_blocks', 's_c_conv', 's_n_conv']:
        sub_df[c] = [ord(s) -65 for s in sub_df[c]]
    c_dict = {s:c for s,c in zip(sub_df['s_blocks'], sub_df['s_c_conv'])}
    n_dict = {s:n for s,n in zip(sub_df['s_blocks'], sub_df['s_n_conv'])}
    dict_dict = {'c':c_dict, 'n':n_dict}
    df['subblock'] = [''.join([code[dict_dict[name[0]][k]]
                               for k in range(len(dict_dict[name[0]]))])
                      for code, name in zip(df['code'], df['name'])]

    # generate sequence as string
    a_and_c = os.path.join(dir, 'expanded_alignment_and_contacts.pkl')
    sample_space, contacts = pickle.load(open(a_and_c, 'rb'))
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
        elif name[0] in set(['c', 'C']):
            seqs.append(''.join(chimera_tools.make_sequence(name_dict[name], c_assignments_dict, sample_space)))
    df['sequence'] = seqs

    # make columns for log properties
    df['log_mKate'] = np.log(df['mKate_mean'])
    df['log_GFP'] = np.log(df['GFP_mean'])
    try:
        df['log_intensity_ratio'] = np.log(df['intensity_ratio'])
        df['log_sum_ratio'] = np.log(df['sum_ratio'])
    except KeyError:
        try:
            df['log_cell_ratio'] = np.log(df['cell_ratio'])
        except:
            pass

    # make binary splits
    split_me = ['mKate_mean', 'GFP_mean', 'sum_ratio',
                'intensity_ratio', 'cell_ratio']
    hows = ['median', 'parents']
    for sp in split_me:
        for how in hows:
            split(df, sp, how=how)

    save_me = os.path.join(dir, args.name.split('/')[0])
    if args.min_cells is not None:
        save_me += '/' + str(args.min_cells) + 'props'
    else:
        save_me += '/props'
    # pickle
    with open(save_me + '.pkl', 'wb') as f:
        pickle.dump(df, f)
    # save as a csv
    with open(save_me + '.csv', 'wb') as f:
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
