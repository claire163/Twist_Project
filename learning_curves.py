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
from sklearn import metrics

def main():
    dir = os.path.dirname(__file__)

    # load the GPModel
    ms = ['2016-05-11/models/2016-06-13__mKate_above_parent_structure.pkl']
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
        ns = range(max_n-6, max_n)
        ns = range(60, max_n)
        ns = ns + [max_n-1]
        #ns = range(10, len(X) - 1, 15)
        for n in ns:
            predicted, actual, R = gptools.cv(X, Y, model, n,
                                              replicates=1000,
                                              keep_inds = ssw)
            Rs.append(R)
            print n, R
#             if model.regr:
#                 Rs.append(np.corrcoef(predicted, actual)[0,1])
#                 print n, Rs[-1]
#             else:
#                 fpr, tpr, _ = metrics.roc_curve(actual, predicted)
#                 auc = metrics.auc(fpr,tpr)
#                 print n, auc

    print Rs
    print ns
    plt.plot(ns, Rs, '.-')
    plt.xlabel('Number maximally informative in training set')
    if model.regr:
        plt.ylabel('R')
    else:
        plt.ylabel('auc')
#     plt.legend(["linear",
#                 'exponential',
#                r"Matern, $\nu=\frac{3}{2}$",
#                 r'Matern, $\nu=\frac{5}{2}$'],
#                loc='best')
    plt.margins(0.02)
    plt.savefig('2016-05-11/plots/mKate_above_parent_lc_avg_zoomed.pdf')
    plt.show()


if __name__ == "__main__":
    main()
