"""
Loads a pickled GPModel and uses it to make predictions
"""
import sys
sys.path.append ('/Users/seinchin/Documents/Caltech/Arnold Lab/Programming tools/GPModel')
import argparse, gpmodel, gpkernel, os, pickle, gptools
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

#     vs = np.linspace(1e-1,5.,100)
#     param = 0
#     res = [log_ML(model, (param,v)) for v in vs]
#     fit = [r[0] for r in res]
#     pen = [r[1] for r in res]
#     plt.plot (vs, fit)
#     plt.xlabel(r'$\sigma_{p}^{2}$')
#     #plt.legend (['fit', 'complexity'])
#     plt.title (r'$\sigma_n^2 = %.5f$' %param)
#     plt.show()

    #for vp in np.linspace(10, 100, 100):
    #    model.logistic_log_ML(vp)
    seq = model.X_seqs.loc[['c4','c36']]
    print model.f_hat
    res = model.predicts(seq)
    print res
#     with open(os.path.join(dir,'29_Sept_pis.txt'),'w') as f:
#         for p in pis:
#             f.write(str(p) + '\n')


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
