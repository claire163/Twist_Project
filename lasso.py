"""
Run Lasso or Bayesian Ridge on a property.

If -a given, use Lasso to select features, then report LOO results
for Bayesian Ridge.

Otherwise, use LOO CV to find best alpha, then use Lasso to select
features and report LOO results for Bayesian Ridge.
"""
# all the useful packages
import numpy as np
import pandas as pd
from sys import exit
from sklearn import linear_model, cross_validation
import sys
sys.path.append('/Users/kevinyang/Documents/Projects/GPModel')
sys.path.append ('/Users/seinchin/Documents/Caltech/Arnold Lab/Programming tools/GPModel')
import chimera_tools, gpkernel, gpmodel
import cPickle as pickle
import matplotlib.pyplot as plt
import os
import argparse
from scipy import stats
import itertools

def LOO(clf, Xs, Ys):
    y = []
    for train, test in cross_validation.LeaveOneOut(len(Xs)):
        X_train = Xs.iloc[train]
        y_train = Ys.iloc[train]
        X_test = Xs.iloc[test]
        clf.fit(X_train,y_train)
        y.append(clf.predict(X_test)[0])
    R = np.corrcoef(np.array(y), Ys)[0][1]
    return R

def Bayesian_Ridge_results(Xs, Ys, subblocks, plot, collapse, y_c):
    clf = linear_model.BayesianRidge()
    if subblocks:
        X_terms = [chimera_tools.get_terms(x) for x in df[inds]['subblock']]
        Xs = chimera_tools.X_from_terms(X_terms, terms)
        print X_terms
    else:
        a_and_c = 'alignment_and_contacts.pkl'
        sample_space, contacts = pickle.load(open(a_and_c))
        Xs, terms = chimera_tools.make_X(df[inds]['sequence'],
                                         sample_space, contacts,
                                         terms=terms,
                                         collapse=args.collapse)
    Xs = pd.DataFrame(Xs, index=Ys.index)
    R = LOO(clf, Xs, Ys)
    print 'R =' %R
    if plot:
        plt.plot(Ys, y,'k.')
        plt.margins(0.02)
        plt.xlabel('actual ' + y_c)
        plt.ylabel('predicted ' +  y_c)
    save_me = args.training.split('/')[0] + '/models/' + y_c + '_' + str(alpha)
    if collapse:
        save_me += '_col'
    if subblocks:
        save_me += '_subblocks'
    plt.savefig(save_me + '_ridge.pdf')
    plt.show()
    with open(save_me + '_lasso_weights.csv', 'w') as f:
        weights.to_csv(f)
    weights = pd.DataFrame()
    weights['weight'] = clf.coef_
    weights['term'] = terms
    with open(save_me + '_ridge_weights.csv', 'w') as f:
        weights.to_csv(f)

def GP_results(kernel_name, Xs, Ys, subblocks, plot,
               write, training, y_c, alpha, terms):
    if kernel_name == 'lin':
        kernel = gpkernel.LinearKernel()
    elif kernel_name == 'SE':
        kernel = gpkernel.SEKernel()
    elif kernel_name == '52':
        kernel = gpkernel.MaternKernel('5/2')
    elif kernel_name == '32':
        kernel = gpkernel.MaternKernel('3/2')
    model = gpmodel.GPModel(kernel)
    if isinstance(Xs, pd.DataFrame):
        Xs.index = Ys.index
    else:
        Xs = pd.DataFrame(Xs, index=Ys.index)
    model.fit(Xs, Ys)
    print model.hypers
    LOOs = model.LOO_res (model.hypers, add_mean=True)
    actual = Ys
    predicted = LOOs['mu']
    var = LOOs['v']
    print 'log_ML = %f' %-model.ML
    print 'log_LOO_P = %f' %-model.log_p
    r1 = stats.rankdata(actual)
    r2 = stats.rankdata(predicted)
    tau = stats.kendalltau(r1, r2).correlation
    R = np.corrcoef(actual, predicted)[0,1]
    print 'tau = %.4f' %tau
    print 'R = %.4f' %R
    print 'Saving model results...'
    if subblocks:
        kernel_name += '_subblocks'
    save_me = [kernel_name, y_c, str(-model.ML),
               str(-model.log_p), str(R), str(tau),
               '', ' '.join(str(model.hypers).split(',')), 'no',
               str(write), str(alpha)]

    with open(training.split('/')[0] + '/models.csv', 'r') as f:
        lines = [','.join(ell.split(',')[0:-2]) for ell in f.readlines()]
    if not ','.join(save_me[0:-2]) in lines:
        with open(training.split('/')[0] + '/models.csv', 'a') as f:
            f.write('\n' + ','.join(save_me))
    if plot:
        plt.plot(Ys, predicted,'k.')
        plt.margins(0.02)
        plt.xlabel('actual ' + y_c)
        plt.ylabel('predicted ' +  y_c)
        plt.show()
    if write:
        name = training.split('/')[0] + '/models/' +\
                    kernel_name + '_' + y_c + '_' + str(alpha)
        if subblocks:
            name += '_subblocks'
        with open(name + '_LOO.txt', 'w') as f:
            f.write(','.join(save_me))
            f.write('\nname,mu,var,y\n')
            for i,n in enumerate(Ys.index):
                f.write (str(n)+','+str(predicted[n]))
                f.write(','+str(var[n])+','+str(Ys.loc[n]))
                f.write('\n')
        if plot:
            plt.savefig(name + '_LOO.pdf')
        model.dump(name + '.pkl')
        with open(name + '_terms.pkl', 'wb') as f:
            pickle.dump(terms, f)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--alpha', type=float, required=False)
    parser.add_argument('-t', '--training',required=True)
    parser.add_argument('-y', '--y_column', required=True)
    parser.add_argument('-c', '--collapse', action='store_true')
    parser.add_argument('-s', '--subblocks', action='store_true')
    parser.add_argument('-k', '--kernel', required=False)
    parser.add_argument('-p', '--plot', action='store_true')
    parser.add_argument('-w', '--write', action='store_true')
    args = parser.parse_args()

    with open(args.training) as f:
        df = pickle.load(f)
    y_c = args.y_column
    inds = ~pd.isnull(df[y_c])
    Ys = df[inds][y_c]
    Ys.index = df[inds]['name']
    X_name = args.training.split('/')[0] + '/X_' + y_c
    if args.collapse:
        X_name += '_col'
    if args.subblocks:
        X_name += '_subblocks'
    print X_name
    print 'Trying to load X...'
    try:
        with open(X_name + '.pkl') as f:
            Xs, terms = pickle.load(f)
            print '\t success!'
    except:
        print 'Building X ...'
        if args.subblocks:
            n_subblocks = len(df.iloc[0]['subblock'])
            sample_space = [('0','1','2') for _ in range(n_subblocks)]
            contacts = [c for c in
                        itertools.combinations(range(n_subblocks), 2)]
            Xs, terms = chimera_tools.make_X(df[inds]['subblock'],
                                             sample_space, contacts)
        else:
            a_and_c = 'alignment_and_contacts.pkl'
            sample_space, contacts = pickle.load(open(a_and_c))
            Xs, terms = chimera_tools.make_X(df[inds]['sequence'],
                                             sample_space, contacts,
                                             collapse=args.collapse)
        Xs = pd.DataFrame(Xs, index=Ys.index)
        with open(X_name + '.pkl', 'wb') as f:
            pickle.dump((Xs, terms), f)

    if args.alpha is None:
        print 'Finding alpha by cross-validation...'
        alphas = np.linspace(0.008, 0.028, num=41)
        Rs = [LOO(linear_model.Lasso(alpha=alpha), Xs, Ys,) for alpha in alphas]
        res = pd.DataFrame()
        res['alpha'] = alphas
        res['R'] = Rs
        print res
        res = res.sort_values('R', ascending=False)
        alpha = res.iloc[0]['alpha']
    else:
        alpha = args.alpha
    print alpha
    if alpha != 0:
        print 'Using Lasso for feature selection...'
        clf = linear_model.Lasso(alpha=alpha)
        clf.fit(Xs, Ys)
        weights = pd.DataFrame()
        weights['weight'] = clf.coef_
        weights['term'] = terms
        weights = weights[~np.isclose(weights['weight'], 0.0)]
        terms = list(weights['term'].values)
        if args.subblocks:
            X_terms = [chimera_tools.get_terms(x) for x in df[inds]['subblock']]
            Xs = chimera_tools.X_from_terms(X_terms, terms)
        else:
            a_and_c = 'alignment_and_contacts.pkl'
            sample_space, contacts = pickle.load(open(a_and_c))
            Xs, terms = chimera_tools.make_X(df[inds]['sequence'],
                                             sample_space, contacts,
                                             terms=terms,
                                             collapse=args.collapse)
    if args.kernel is not None:
        GP_results(args.kernel, Xs, Ys, args.subblocks,
                   args.plot, args.write, args.training, y_c, alpha, terms)
    else:
        Bayesian_Ridge_results(Xs, Ys, args.subblocks,
                               args.plot, args.collapse, y_c)

if __name__ == '__main__':
    main()




