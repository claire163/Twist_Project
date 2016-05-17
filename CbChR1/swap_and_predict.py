import sys
sys.path.append('/Users/kevinyang/Documents/Projects/GPModel')
sys.path.append ('/Users/seinchin/Documents/Caltech/Arnold Lab/Programming tools/GPModel')
sys.path.append ('/Users/seinchin/Documents/Caltech/Arnold Lab/Programming tools/Twist Project')

import argparse, gpmodel, gpkernel, os, gptools
import cPickle as pickle
from chimera_tools import substitute_blocks, load_assignments
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import datetime
import itertools

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--model',required=True)
    parser.add_argument('-p', '--penalty', default=0.0, type=float)
    parser.add_argument('-w', '--write', required=False)

    args = parser.parse_args()

    print 'Loading original sequences...'
    df = pd.read_csv('CbChR1_constructs.txt')
    print 'Making single block swaps'
    # load the assignment dicts for contiguous and non-contiguous
    c_assignments_file = '../clibrary.output'
    n_assignments_file = '../nlibrary.output'
    c_assignments_dict = load_assignments(c_assignments_file)
    n_assignments_dict = load_assignments(n_assignments_file)
    a_and_c = '../alignment_and_contacts.pkl'
    sample_space, contacts = pickle.load(open(a_and_c))
    n_parents = len(sample_space[0])
    n_blocks = max(c_assignments_dict.values())+1
    names = []
    sequences = []
    mutations = []
    # make single swaps
    for i in range(n_blocks):
        for j in range(n_parents):
            for k in range(2):
                st = str(n_parents+k)
                name = ''.join(i * [st] + [str(j)] + max(0, n_blocks-i-1) * [st])
                names.append('c' + name)
                names.append('n' + name)
                if k == 0:
                    ori = 'CbChR1'
                else:
                    ori = 'CsCbChR1'
                original = df[df['name'] ==ori]['sequence'].values[0]
                for a_d in [c_assignments_dict, n_assignments_dict]:
                    sequences.append(substitute_blocks(original,
                                                       [(j,i)],
                                                       a_d,
                                                       sample_space))
                    m = sum([1 for o,s in zip(original, sequences[-1])
                             if o != s])
                    mutations.append(m)
    df = df.append(pd.DataFrame({'name':names,
                                 'm':mutations,
                                 'sequence':sequences}))
    X = [[s for s in sequence] for sequence in df['sequence']]
    X = pd.DataFrame(X, index = df['name'])
    print 'Loading the model...'
    model = gpmodel.GPModel.load(args.model)
    print 'Making predictions'
    preds = np.array(model.predicts(X))
    df['mean'] = preds[:,0]
    df['variance'] = preds[:,1]
    df['UB2'] = df['mean'] + 2 * np.sqrt(df['variance'])
    df = df.sort_values('UB2', ascending=False)
    df = df.sort_values('m')
    pareto = []
    best = df.iloc[0]['UB2']
    for row in df.itertuples(index=False, name='row'):
        if row.UB2 > best:
            best = row.UB2
            pareto.append(True)
        else:
            pareto.append(False)
    df['pareto'] = pareto
    print df[df['pareto']][['name', 'UB2', 'm']]
    df = df.sort_values('pareto', ascending=False)
    df = df[['name', 'UB2', 'm', 'pareto', 'mean', 'variance', 'sequence']]
    if args.write:
        print 'Writing results to file...'
        df.to_csv(args.write, index=False)


if __name__=="__main__":
    main()
