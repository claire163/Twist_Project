# all the useful packages
import numpy as np
import pandas as pd
from sys import exit
from sklearn import cross_validation
from sklearn import grid_search
from sklearn import ensemble
from sklearn import linear_model
import sys
sys.path.append('/Users/kevinyang/Documents/Projects/GPModel')
sys.path.append ('/Users/seinchin/Documents/Caltech/Arnold Lab/Programming tools/GPModel')
import chimera_tools
import cPickle as pickle
import matplotlib.pyplot as plt
import os
import datetime
import argparse

def normalize(data):
    """
        Normalizes the elements in data by subtracting the mean and dividing
        by the standard deviation.

        Parameters:
            data (pd.Series)

        Returns:
           normed
        """
    m = data.mean()
    s = data.std()
    return (data-m) / s



parser = argparse.ArgumentParser()
parser.add_argument('-a', '--params', required=False, nargs='*')
parser.add_argument('-m', '--method', required=True)
parser.add_argument('-t', '--training',required=True)
parser.add_argument('-L', '--LOO', action='store_true')
parser.add_argument('-y', '--y_column', required=True)
parser.add_argument('-p', '--plot', action='store_true')
parser.add_argument('-d', '--drop', required=False, type=float)


args = parser.parse_args()
if args.params is None:
    args.params = {}
else:
    args.params = {p.split(':')[0]:float(p.split(':')[1]) for p in args.params}
with open(args.training + '/props.pkl') as f:
    df = pickle.load(f)
y_c = args.y_column
inds = ~pd.isnull(df[y_c])
Ys = df[inds][y_c]
Ys.index = range(len(Ys))
dt = str(datetime.date.today())
X_name = 'sk_methods/' + args.training + '_' + y_c
funct_dict = {'Lasso':linear_model.Lasso,
              'BayesianRidge':linear_model.BayesianRidge}

print 'Trying to load X...'
try:
    with open(X_name + '.pkl') as f:
        Xs, terms = pickle.load(f)
        print '\t success!'
except:
    print 'Building X ...'
    a_and_c = 'alignment_and_contacts.pkl'
    sample_space, contacts = pickle.load(open(a_and_c))
    Xs, terms = chimera_tools.make_X(df[inds]['sequence'], sample_space, contacts)
    Xs = pd.DataFrame(Xs, index=Ys.index)
    with open(X_name + '.pkl', 'wb') as f:
        pickle.dump((Xs, terms), f)

clf = funct_dict[args.method](**args.params)
if args.LOO:
    y = []
    for train, test in cross_validation.LeaveOneOut(len(Xs)):
        X_train = Xs.loc[train]
        y_train = Ys.loc[train]
        X_test = Xs.loc[test]
        clf.fit(X_train,y_train)
        y.append(clf.predict(X_test)[0])
    print 'R =', np.corrcoef(np.array(y), Ys)[0][1]
    plt.plot(Ys, y,'k.')
    plt.margins(0.02)
    plt.xlabel('actual ' + y_c)
    plt.ylabel('predicted ' +  y_c)
    #plt.savefig(X_name + '_' + args.method + '.pdf')
    plt.show()
else:
    clf.fit(Xs, Ys)
    weights = pd.DataFrame()
    weights['weight'] = clf.coef_
    weights['term'] = terms
    weights = weights[~np.isclose(weights['weight'], 0.0)]
    print weights
    with open(X_name + '_' + args.method + '_weights.csv', 'w') as f:
        weights.to_csv(f)


# collapse redundant weights
# use these Xs with GP

