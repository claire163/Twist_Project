'''
Reads in a model and a validation set, makes predictions, and compares the results.

'''
import sys
sys.path.append('/Users/kevinyang/Documents/Projects/GPModel')
sys.path.append ('/Users/seinchin/Documents/Caltech/Arnold Lab/Programming tools/GPModel')
import argparse, gpmodel, gpkernel, os, gptools
import dill as pickle
from chimera_tools import *
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import datetime
from scipy import stats
from collections import namedtuple
from sklearn import metrics
import seaborn as sns


def main ():
    dir = os.path.dirname(__file__)
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--model',required=True)
    parser.add_argument('-v', '--validate',required=True)
    parser.add_argument('-y', '--y_column',required=True)
    parser.add_argument('-p', '--plot', action='store_true')
    parser.add_argument('-w', '--out_file', required=False)


    args = parser.parse_args()

    print 'Loading model...'
    model = gpmodel.GPModel.load(args.model)
    try:
        with open(args.model.split('.pkl')[0] + '_terms.pkl', 'rb') as f:
            terms = pickle.load(f)
            has_terms = True
    except IOError:
        has_terms = False

    print 'Loading validation set...'
    with open (os.path.join(dir,args.validate),'r') as f:
        df = pickle.load(f)
    parents = ['cschrimson', 'c1c2', 'cheriff']
    df = df[~df['name'].isin(parents)]

    print 'Making predictions...'
    X_seqs = [list(seq) for seq in df['sequence']]
    X_seqs = pd.DataFrame(X_seqs, index=df['name'])
    if has_terms:
        a_and_c = 'alignment_and_contacts.pkl'
        sample_space, contacts = pickle.load(open(a_and_c))
        X, terms = make_X(df['sequence'], sample_space, contacts,
                          terms=terms, collapse=False)
        X = pd.DataFrame(X, index=df['name'])
        data = pd.DataFrame(model.predicts(X),
                            index=df['name'])
    else:
        data = pd.DataFrame(model.predicts(X_seqs),
                            index=X_seqs.index)
    print 'Computing metrics...'
    if model.regr:
        data.columns = ['mu', 'variance']
        y = df[args.y_column]
        y.index = data.index
        data['y'] = y
        print data
        r1 = stats.rankdata(data['y'])
        r2 = stats.rankdata(data['mu'])
        tau = stats.kendalltau(r1, r2).correlation
        R = np.corrcoef(data['y'], data['mu'])[0,1]
        print 'tau = %.4f' %tau
        print 'R = %.4f' %R
    else:
        data.columns = ['pi', 'f_bar', 'variance']
        y = df[args.y_column]
        y.index = data.index
        data['y'] = y
        print data
        fpr, tpr, _ = metrics.roc_curve(data['y'], data['pi'])
        auc = metrics.auc(fpr,tpr)
        print 'AUC = %.4f' %auc
    # plot
    # write results
    if args.out_file is not None:
        print 'Writing results...'
        data.to_csv(args.out_file)
        with open(args.out_file, 'a') as f:
            if model.regr:
                f.write('# tau = %.4f\n' %tau)
                f.write('# R = %.4f' %R)
            else:
                f.write('# AUC = %.4f' %auc)


    if args.plot:
        print 'Plotting...'
        if model.regr:
            std = np.sqrt(data['variance'])
            plt.plot(data['y'], data['mu'], 'o')
            plt.margins(0.02)
            plt.xlabel('actual ' + args.y_column)
            plt.ylabel('predicted' + args.y_column)
        else:
            auc = gptools.plot_ROC(data['y'], data['pi'])
        if args.out_file is not None:
            plt.savefig(args.out_file.split('.txt')[0] + '.pdf')
        else:
            plt.show()



if __name__ == "__main__":
    main()
