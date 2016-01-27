# all the useful packages
import numpy as np
import pandas as pd
from sys import exit
from sklearn import cross_validation
from sklearn import ensemble
import chimera_tools
import cPickle as pickle
import matplotlib.pyplot as plt
import os

def code_to_x (code):
    x = np.zeros(len(code)*3)
    ind = 0
    x[0]=1
    for c in code:
        ind += int(c) - 1
        x[ind] = 1
        ind += 4 - int(c)
    return x

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


with open('2015-12-07_all_res.pkl') as f:
    df = pickle.load(f)

y_c = 'log_gfp'
a_and_c = 'alignment_and_contacts.pkl'
sample_space, contacts = pickle.load(open(a_and_c))

inds = ~np.isnan(df[y_c])
terms = chimera_tools.contacting_terms(sample_space, contacts)
Xs = chimera_tools.make_contact_X(df[inds]['sequence'], contacts, terms)
Ys = df[inds][y_c]

X_train, X_test, y_train, y_test = \
    cross_validation.train_test_split(Xs, Ys, test_size=0.2)


clf = ensemble.RandomForestRegressor(n_estimators=200)
clf.fit(X_train,y_train)
y = clf.predict(X_test)



print 'R =', np.corrcoef(y, y_test)[0][1]
print 'validation score:', clf.score(X_test, y_test)
plt.plot(y_test, y,'k.')
# plt.plot(res[y_c], res['predicted'], 'k.', alpha=0.3)
# plt.xlabel('actual ' + y_c)
# plt.ylabel('predicted ' +  y_c)
# #plt.savefig('plots/2015-12-07_contig_lasso_validation.pdf')
plt.show()

