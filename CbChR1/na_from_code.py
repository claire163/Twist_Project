import sys
sys.path.append('/Users/kevinyang/Documents/Projects/GPModel')
sys.path.append ('/Users/seinchin/Documents/Caltech/Arnold Lab/Programming tools/GPModel')
sys.path.append ('/Users/seinchin/Documents/Caltech/Arnold Lab/Programming tools/Twist Project')

import argparse
import cPickle as pickle
from chimera_tools import make_sequence, load_assignments, translate
import pandas as pd

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--codes', required=True, nargs='*')
    parser.add_argument('-w', '--out_file', required=True)
    args = parser.parse_args()

    print 'Loading sample spaces and assignments...'
    # load the assignment dicts for contiguous and non-contiguous
    c_assignments_file = '../clibrary.output'
    n_assignments_file = '../nlibrary.output'
    c_assignments_dict = load_assignments(c_assignments_file)
    n_assignments_dict = load_assignments(n_assignments_file)
    assign = {'c':c_assignments_dict, 'n': n_assignments_dict}
    a_and_c = '../alignment_and_contacts.pkl'
    na_ss_file = '../na_sample_space.pkl'
    aa_ss, _ = pickle.load(open(a_and_c))
    n_parents = len(aa_ss[0])
    n_blocks = max(c_assignments_dict.values()) + 1
    sequences = []
    with open('../parents_na.txt', 'r') as f:
        df = pd.read_csv(f)
        df['sequence'] = [seq.lower() for seq in df['sequence']]
    try:
        na_ss = pickle.load(open(na_ss_file))
    except IOError:
        # chop into codons
        na_ss = [tuple([row['sequence'][j:j+3] for i, row in df.iterrows()])
                 for j in range(0,len(df.iloc[0,1]),3)]
        new_na_ss = []
        # make conserved residues use c1c2 encoding
        for aas, nas in zip(aa_ss, na_ss):
            if aas[0] == aas[1] and aas[0] == aas[2]:
                new_na_ss.append((nas[1], nas[1], nas[1], nas[3], nas[4]))
            else:
                new_na_ss.append(nas)
        na_ss = new_na_ss
        with open(na_ss_file, 'wb') as f:
            pickle.dump(na_ss, f)
    assert ''.join([aa[3] for aa in na_ss]) == df[df['name']
                                                  =='cbchr1']['sequence'].values[0]
    assert ''.join([aa[4] for aa in na_ss]) == df[df['name']
                                                  =='cscbchr1']['sequence'].values[0]
    print 'Generating sequences...'
    new_df = pd.DataFrame(args.codes, columns=['code'])
    new_df['na_seq'] = [''.join(make_sequence(code[1::], assign[code[0]], na_ss))
                        for code in args.codes]
    new_df['aa_seq'] = [translate(na) for na in new_df['na_seq']]
    new_df.to_csv(args.out_file, index=False)

if __name__=="__main__":
    main()
