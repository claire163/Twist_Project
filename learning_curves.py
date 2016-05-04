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
    ms = ['2016-02-17_models/2016-03-25__log_mKate_structure_dict.pkl',
         '2016-02-17_models/2016-03-25__log_mKate_SEStructure_dict.pkl',
         '2016-02-17_models/2016-03-25__log_mKate_32Structure_dict.pkl',
         '2016-02-17_models/2016-03-25__log_mKate_52Structure_dict.pkl']
         #'2016-02-02__log_mKate_SEStructure_dict.pkl']
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
        mi = list(set(X.index) - set(ssw))
        X_ssw = X.loc[ssw]
        X_mi = X.loc[mi]
        Y_ssw = Y.loc[ssw]
        Y_mi = Y.loc[mi]
        Rs = []
        ns = range(0, len(mi), 15)
        ns = ns + [53]
        for n in ns[::-1]:
            print n
            predicted, actual = gptools.cv(X_mi, Y_mi, model, n,
                                           replicates=100,
                                          X_always=X_ssw,
                                          Y_always=Y_ssw)
            Rs.append(np.corrcoef(predicted, actual)[0,1])
        plt.plot(ns, Rs[::-1], '.-')
    plt.xlabel('Number maximally informative in training set')
    plt.ylabel('R')
    plt.legend(["linear",
                'exponential',
               r"Matern, $\nu=\frac{3}{2}$",
                r'Matern, $\nu=\frac{5}{2}$'],
               loc='best')
    plt.margins(0.02)
    plt.savefig('plots/log_mKate_mean_learning_curves_3.pdf')
    plt.show()


if __name__ == "__main__":
    main()
