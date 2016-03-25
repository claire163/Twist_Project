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
import datetime

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

y_c = 'log_kinetics'
dt = str(datetime.date.today())

n = 100
name = dt + '_struct_X_Y_' + y_c
print "Trying to load X and Y..."
try:
    with open(name + '.pkl') as f:
        Xs, Ys = pickle.load(f)
        print '\t success!'
except:
    print 'Building X and Y...'
    with open('2016-03-21_data/props.pkl') as f:
        df = pickle.load(f)
    a_and_c = 'alignment_and_contacts.pkl'
    sample_space, contacts = pickle.load(open(a_and_c))

    inds = ~np.isnan(df[y_c])
    terms = chimera_tools.contacting_terms(sample_space, contacts)
    Xs = chimera_tools.make_contact_X(df[inds]['sequence'], contacts, terms)
    Ys = df[inds][y_c]
    Ys.index = range(len(Ys))
    Xs = pd.DataFrame(Xs, index=Ys.index)
    with open(name + '.pkl', 'wb') as f:
        pickle.dump((Xs, Ys), f)
clf = ensemble.RandomForestRegressor(n_estimators=n)

y = []
print 'Generating LOO results...'
for train, test in cross_validation.LeaveOneOut(len(Xs)):
    X_train = Xs.loc[train]
    y_train = Ys.loc[train]
    X_test = Xs.loc[test]
    clf.fit(X_train,y_train)
    y.append(clf.predict(X_test)[0])

with open('LOO_results/'+dt+'_random_forest_'+str(n)+'_'+y_c+'_LOO.txt','w') as f:
    f.write('predicted,actual\n')
    for p,a in zip(y, Ys):
        f.write (str(p) + ',' + str(a))
        f.write('\n')

print 'R =', np.corrcoef(np.array(y), Ys)[0][1]
plt.plot(Ys, y,'k.')
plt.margins(0.02)
plt.xlabel('actual ' + y_c)
plt.ylabel('predicted ' +  y_c)
plt.savefig('plots/' + dt + '_random_forest_' + str(n) + '_' + y_c + '.pdf')
plt.show()

