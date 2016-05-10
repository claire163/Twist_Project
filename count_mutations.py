""" Count number of mutations in each chimera read in."""

import sys, os
import dill as pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def count_changes(str_1, str_2):
    return sum([1 for a,b in zip(str_1, str_2) if a!=b])

dir = os.path.dirname(__file__)

def main():
    try:
        n_mutations = pd.read_csv('chimera_mutations.txt')
    except:
        sequences = pd.read_csv('all_chimeras.txt')
        sequences.index = xrange(len(sequences))
        codes = sequences.iloc[:,0]
        sequences = sequences.drop(sequences.columns[[0]], axis=1)
        sequences = pd.DataFrame([''.join(s) for _,s in sequences.iterrows()],
                                 index=codes.values,
                                columns=['sequence'])
        parent_inds = ['n0000000000', 'n1111111111', 'n2222222222']
        n_mutations = [min([count_changes(sequences.loc[p,'sequence'], s)
                            for p in parent_inds])
                       for _,s in sequences['sequence'].iteritems()]
        n_mutations = pd.DataFrame(n_mutations,
                                  columns=['n_mutations'])
        n_mutations['code'] = codes
        n_mutations = n_mutations[['code', 'n_mutations']]
        n_mutations.to_csv('chimera_mutations.txt', index=False)
    print np.mean(n_mutations['n_mutations'])
    print np.median(n_mutations['n_mutations'])
    df = pd.read_csv('2016-05-04_predictions/GFP_above_parent_SEStructure.txt')
    plt.plot(n_mutations['n_mutations'], df['pi'], '.')
    plt.show()





if __name__=="__main__":
    main()
