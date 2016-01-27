"""
Loads a pickled GPModel and uses it to make predictions (or whatever else)
"""
import sys
sys.path.append ('/Users/seinchin/Documents/Caltech/Arnold Lab/Programming tools/GPModel')
import argparse, gpmodel, gpkernel, os, gptools
import cPickle as pickle
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math
from scipy import stats
from scipy.optimize import minimize

def main():
    dir = os.path.dirname(__file__)
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--model',required=True)

    # load the GPModel
    args = parser.parse_args()
    ms = ['2016-01-06__log_mKate_structure_dict.pkl',
         '2016-01-06__log_mKate_SEStructure_dict.pkl']
    for m in ms:
    print 'loading model...'
    model = gpmodel.GPModel.load(ms)
    print 'doing cv ...'
    X = model.X_seqs
    Y = model.Y
    taus = []
    Rs = []
    ns = [50, 75, 100, 117]
    for n in ns:
        print n
        predicted, actual = gptools.cv(X, Y, model, n, 1000)
        r1 = stats.rankdata(actual)
        r2 = stats.rankdata(predicted)
        taus.append(stats.kendalltau(r1, r2).correlation)
        Rs.append(np.corrcoef(predicted, actual)[0,1])
    print taus, Rs
    plt.plot(ns, taus, '.-', ns, Rs, '.-')
    plt.xlabel('training set size')
    plt.legend(["structure Kendall's Tau",
                'structure R',
               "SEStructure Kendall's Tau",
                'SEStructure R'],
               loc='best')
    plt.ylim([0,1])
    plt.savefig('2016-01-06_log_mKate_mean_structure_learning_curve.pdf')
    plt.show()


# helper function for plotting fit and penalty terms of ML
def log_ML (model,variances):
        """ Returns the negative log marginal likelihood.

        Parameters:
            variances (iterable): var_n and var_p

        Uses RW Equation 5.8
        """
        Y_mat = np.matrix(model.Y)
        var_n,var_p = variances
        K_mat = np.matrix (model.K)
        Ky = K_mat*var_p+np.identity(len(K_mat))*var_n
        L = np.linalg.cholesky (Ky)
        alpha = np.linalg.lstsq(L.T,np.linalg.lstsq (L, np.matrix(Y_mat).T)[0])[0]
        fit =  (0.5*Y_mat*alpha).item()
        penalty = sum([math.log(l) for l in np.diag(L)])

        return (fit, penalty)

if __name__ == "__main__":
    main()
