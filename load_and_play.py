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
from scipy.optimize import minimize

def main():
    dir = os.path.dirname(__file__)
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--model',required=True)

    # load the GPModel
    args = parser.parse_args()
    with open(os.path.join(dir,args.model),'r') as m_file:
        model = pickle.load(m_file)
    print model.var_n, model.var_p, model.ML

    print model.K.min().min()

#     model.normed_Y = model.Y
#     print minimize(model.log_ML,
#                    (1.0, 1.0),
#                    method='L-BFGS-B',
#                    bounds=[(1e-5,None),(1e-5,None)])
#     res = minimize(model.LOO_MSE, (0.74, 0.4),
#                   method='L-BFGS-B',
#                   bounds=[(1e-5,None),(1e-5,None)])
#     print res

#     print model.LOO_log_p(res['x'])
#     print model.log_ML(res['x'])
#     print minimize(model.log_ML, (0.74, 0.4),
#                   method='L-BFGS-B',
#                   bounds=[(1e-5,None),(1e-5,None)])#     seq = model.X_seqs.loc[['n72', 'n47']]
#     res = model.predicts(seq)
#     print model.Y['n72'], model.Y['n47']
#     print res

#     print model.var_p, model.var_n
#     print model.ML
#     gptools.plot_ML_parts (model, (3.0,20.0),
#                           n=10)
#     plt.show()

#     vars = (10.0)

#     print minimize(model.logistic_log_ML,
#                    vars,
#                    bounds=[(1e-5,None)],
#                    method='L-BFGS-B'
#                   )
#     exit('')
#     print minimize(m.logistic_log_ML,
#                                     10.,
#                                     bounds=[(1e-4, None)])
#     res = plot_ML_parts(m, ranges=([10.0],[0.001,10.0]),
#                         n=100, lab='Normed no zero mKate_mean',
#                        plots=['complexity'])

#     res = plot_ML_contour(m, ranges=([0.72,0.76],[0.0002,.0003]),
#                         n=100, lab='Normed mKate_mean',
#                          n_levels=20)
#     save_as = '2015-11-6_ML_plot_4.pdf'
#     plt.savefig(save_as)
#     plt.show()
#     vars = (2.833, 0.1965)
#     print log_marginal_likelihood(m,vars)
    #print log_ML(m,vars)


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
