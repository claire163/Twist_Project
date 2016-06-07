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
from sklearn import linear_model, cross_validation
import sys
sys.path.append('/Users/kevinyang/Documents/Projects/GPModel')
sys.path.append ('/Users/seinchin/Documents/Caltech/Arnold Lab/Programming tools/GPModel')
import chimera_tools, gpmodel, gpkernel
import cPickle as pickle
import matplotlib.pyplot as plt
import argparse
from scipy import stats

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
parser.add_argument('-a', '--alpha', type=float, required=False)
parser.add_argument('-t', '--training',required=True)
parser.add_argument('-y', '--y_column', required=True)
parser.add_argument('-k', '--kernel', required=True)
parser.add_argument('-p', '--plot', action='store_true')
parser.add_argument('-w', '--write', action='store_true')
args = parser.parse_args()


with open(args.training) as f:
    df = pickle.load(f)
y_c = args.y_column
inds = ~pd.isnull(df[y_c])
Ys = df[inds][y_c]
Ys.index = df[inds]['name']
X_name = args.training.split('/')[0] + '/' + 'X_' + y_c

print 'Trying to load X...'
try:
    with open(X_name + '.pkl') as f:
        Xs, terms = pickle.load(f)
        print '\t success!'
except:
    print 'Building X ...'
    a_and_c = 'alignment_and_contacts.pkl'
    sample_space, contacts = pickle.load(open(a_and_c))
    Xs, terms = chimera_tools.make_X(df[inds]['sequence'],
                                     sample_space, contacts,
                                     collapse=False)
    Xs = pd.DataFrame(Xs, index=Ys.index)
    with open(X_name + '.pkl', 'wb') as f:
        pickle.dump((Xs, terms), f)

if args.alpha is None:
    alphas = np.linspace(0.008, 0.023, num=31)
    Rs = []
    for alpha in alphas:
        clf = linear_model.Lasso(alpha=alpha)
        y = []
        for train, test in cross_validation.LeaveOneOut(len(Xs)):
            X_train = Xs.iloc[train]
            y_train = Ys.iloc[train]
            X_test = Xs.iloc[test]
            clf.fit(X_train,y_train)
            y.append(clf.predict(X_test)[0])
        R = np.corrcoef(np.array(y), Ys)[0][1]
        Rs.append(R)
    res = pd.DataFrame()
    res['alpha'] = alphas
    res['R'] = Rs
    print res
    res = res.sort_values('R', ascending=False)
    alpha = res.iloc[0]['alpha']

else:
    alpha = args.alpha
print alpha
clf = linear_model.Lasso(alpha=alpha)
clf.fit(Xs, Ys)
weights = pd.DataFrame()
weights['weight'] = clf.coef_
weights['term'] = terms
weights = weights[~np.isclose(weights['weight'], 0.0)]
terms = list(weights['term'].values)
if args.kernel == 'lin':
    kernel = gpkernel.LinearKernel()
elif args.kernel == 'SE':
    kernel = gpkernel.SEKernel()
elif args.kernel == '52':
    kernel = gpkernel.MaternKernel('5/2')
elif args.kernel == '32':
    kernel = gpkernel.MaternKernel('3/2')
model = gpmodel.GPModel(kernel)
a_and_c = 'alignment_and_contacts.pkl'
sample_space, contacts = pickle.load(open(a_and_c))
Xs, terms = chimera_tools.make_X(df[inds]['sequence'],
                                 sample_space, contacts,
                                 terms=terms,
                                 collapse=False)
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
save_me = [args.kernel, args.y_column, str(-model.ML),
           str(-model.log_p), str(R), str(tau),
           '', ' '.join(str(model.hypers).split(',')), 'no',
           str(args.write), str(alpha)]

with open(args.training.split('/')[0] + '/models.csv', 'r') as f:
    lines = [','.join(ell.split(',')[0:-2]) for ell in f.readlines()]
if not ','.join(save_me[0:-2]) in lines:
    with open(args.training.split('/')[0] + '/models.csv', 'a') as f:
        f.write('\n' + ','.join(save_me))
if args.plot:
    plt.plot(Ys, predicted,'k.')
    plt.margins(0.02)
    plt.xlabel('actual ' + y_c)
    plt.ylabel('predicted ' +  y_c)
if args.write:
    name = args.training.split('/')[0] + '/models/' +\
                args.kernel + '_' + y_c + '_' + str(alpha)
    if args.plot:
        plt.savefig(name + '_LOO.pdf')
        with open(name + '_LOO.txt', 'w') as f:
            f.write(','.join(save_me))
            f.write('\nname,mu,var,y\n')
            for i,n in enumerate(Ys.index):
                f.write (str(n)+','+str(predicted[n]))
                f.write(','+str(var[n])+','+str(Ys.loc[n]))
                f.write('\n')
        plt.show()
    model.dump(name + '.pkl')
    with open(name + '_terms.pkl', 'wb') as f:
        pickle.dump(terms, f)
