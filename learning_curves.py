"""
Loads a pickled GPModel and uses it to make predictions (or whatever else)
"""
import sys
sys.path.append ('/Users/seinchin/Documents/Caltech/Arnold Lab/Programming tools/GPModel')
sys.path.append('/Users/kevinyang/Documents/Projects/GPModel')

import gpmodel, gpkernel, os, gptools
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy import stats
from scipy.optimize import minimize

def main():
    dir = os.path.dirname(__file__)

    # load the GPModel
    ms = ['2016-02-02__log_mKate_structure_dict.pkl']
         #'2016-02-02__log_mKate_SEStructure_dict.pkl']
    for m in ms:
        print 'loading model...'
        model = gpmodel.GPModel.load(m)
        print 'doing cv ...'
        X = model.X_seqs
        Y = model.Y
        taus = []
        Rs = []
        ns = range(10,117)
        for n in ns:
            print n
            predicted, actual = gptools.cv(X, Y, model, n, 10)
            r1 = stats.rankdata(actual)
            r2 = stats.rankdata(predicted)
            taus.append(stats.kendalltau(r1, r2).correlation)
            Rs.append(np.corrcoef(predicted, actual)[0,1])
            #print taus, Rs
        plt.plot(ns, taus, '.-', ns, Rs, '.-')
    plt.xlabel('training set size')
#     plt.legend(["structure Kendall's Tau",
#                 'structure R',
#                "SEStructure Kendall's Tau",
#                 'SEStructure R'],
#                loc='best')
    plt.ylim([0,1])
    plt.savefig('2016-02-2_log_mKate_mean_learning_curves_1.pdf')
    plt.show()


if __name__ == "__main__":
    main()
