"""
Makes a learning curve.
"""
import sys
sys.path.append ('/Users/seinchin/Documents/Caltech/Arnold Lab/Programming tools/GPModel')
sys.path.append('/Users/kevinyang/Documents/Projects/GPModel')

import gpmodel, gpkernel, os, gptools
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy import stats

def main():
    dir = os.path.dirname(__file__)

    # load the GPModel
    ms = ['2016-05-11/models/2016-06-13__log_mKate_structure.pkl']
    for m in ms:
        print 'loading model...'
        model = gpmodel.GPModel.load(m)
        print 'doing cv ...'
        X = model.X_seqs
        Y = model.Y
        ssw = []
        for i in X.index:
            try:
                j = int(i[1::])
                if j < 21:
                    ssw.append(i)
                elif j > 36 and j < 74:
                    ssw.append(i)
            except:
                ssw.append(i)
        Rs = []
        max_n = len(X) - len(ssw)
        ns = range(0, max_n, 15)
        ns = ns + [max_n-1]
        for n in ns[::-1]:
            print n
            predicted, actual = gptools.cv(X, Y, model, n,
                                           replicates=1,
                                          keep_inds = ssw)
            Rs.append(np.corrcoef(predicted, actual)[0,1])
    print Rs
#         plt.plot(ns, Rs[::-1], '.-')
#     plt.xlabel('Number maximally informative in training set')
#     plt.ylabel('R')
#     plt.legend(["linear",
#                 'exponential',
#                r"Matern, $\nu=\frac{3}{2}$",
#                 r'Matern, $\nu=\frac{5}{2}$'],
#                loc='best')
#     plt.margins(0.02)
#     plt.savefig('plots/log_mKate_mean_learning_curves_3.pdf')
#     plt.show()


if __name__ == "__main__":
    main()
