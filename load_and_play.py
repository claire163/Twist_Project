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

def main():
    dir = os.path.dirname(__file__)
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--model',required=True)

    # load the GPModel
    args = parser.parse_args()
    with open(os.path.join(dir,args.model),'r') as m_file:
        model = pickle.load(m_file)


    seq = model.X_seqs.loc[['c5', 'n47']]
    res = model.predicts(seq)
    print res


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
