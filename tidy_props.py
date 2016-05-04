"""
Reads in the measurements, adds columns for the chimera codes and sequence,
and writes it to a csv.
"""

import pandas as pd
import argparse
import os
import chimera_tools
import pickle
from sys import exit
import numpy as np


def main():
    dir = os.path.dirname(__file__)
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--name',required=True)

    args = parser.parse_args()

    df = pd.read_excel(args.name)
    df.columns = [s.lower() for s in df.columns]
    # make names lowercase
    df['name'] = [s.lower() for s in df['name']]

    df['kinetics_no_high'] = df[df.kinetics < 200]['kinetics']

    # calculate logs
    columns = df.columns
    for c in columns:
        try:
            if c[0:4] == 'peak':
                df['log_'+c] = np.log(-df[c])
            else:
                df['log_'+c] = np.log(df[c])
        except:
            pass

    # get the codes from the names
    dict_df = pd.read_excel('dict.xlsx')
    code_dict = {n.lower():c for n,c in zip(dict_df['name'], dict_df['code'])}
    codes = [code_dict[n] for n in df['name']]
    df['code'] = codes

    # zero index codes
    df['code'] = [chimera_tools.zero_index(c) for c in df['code']]

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
            seqs.append(''.join(chimera_tools.make_sequence(name_dict[name],
                                                            n_assignments_dict,
                                                            sample_space)))
        else:
            seqs.append(''.join(chimera_tools.make_sequence(name_dict[name],
                                                            c_assignments_dict,
                                                            sample_space)))
    df['sequence'] = seqs
    # pickle
    with open(os.path.join(dir, args.name.split('/')[0]
                           + '/props.pkl'), 'wb') as f:
        pickle.dump(df, f)
    # save as a csv
    with open(os.path.join(dir, args.name.split('/')[0]
                           + '/props.csv'), 'wb') as f:
        df.to_csv(f)

if __name__ == "__main__":
    main()
