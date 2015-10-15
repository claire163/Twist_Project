""" read in LOO data and plot it"""
import matplotlib.pyplot as plt
import pandas as pd

import sys
sys.path.append ('/Users/seinchin/Documents/Caltech/Arnold Lab/Programming tools/GPModel')
import gpmodel, gpkernel, os, gptools, pickle
import numpy as np
import datetime

dir = os.path.dirname(__file__)
fname = '2015-10-09_test_expression_structure_LOO.txt'
with open (os.path.join(dir,fname),'r') as f:
    df = pd.read_csv(f,header=None)

df.columns = ['name','actual','predicted']

gptools.plot_ROC (df['actual'], df['predicted'], file_name='2015-10-14_exp_struct_kernel_LOO_ROC.pdf',title='Expression structure LOO')
plt.show()



# gfp_name = '2015-10-13__gfp_structure_LOO'
# mKate_name = '2015-10-13_zeroed_mKate_mean_structure_LOO'

# with open (os.path.join(dir,'2015-10-07_all_res.pkl'),'r') as f:
#     train_df = pickle.load(f)

# gfp_df = pd.read_csv(os.path.join(dir,gfp_name+'.txt'),header=None)
# gfp_df.columns = ['name','a_gfp','p_gfp', 'gfp_var']
# gfp_df.index = gfp_df['name']

# mKate_df = pd.read_csv(os.path.join(dir,mKate_name+'.txt'),header=None)
# mKate_df.columns = ['name','a_mKate','p_mKate', 'mKate_var']
# mKate_df.index = mKate_df['name']

# train_inds = ~pd.isnull(train_df['mKate_mean'])
# first = mKate_df['a_mKate'] > 0
# second = mKate_df['p_mKate'] > 0

# predictions = []
# actual = []
# for n in mKate_df[first & second].index:
#     predictions.append(gfp_df.loc[n]['p_gfp']/mKate_df.loc[n]['p_mKate'])
#     a = train_df[train_df['name'] == n]['sum_ratio']
#     a = a[a.index[0]]
#     actual.append(a)


# plt.plot (actual, predictions, '.')
# plt.plot((0,6),(0,6))
# plt.xlabel ('actual sum_ratio')
# plt.ylabel ('predicted sum_ratio')
# plt.savefig('2015-10-13-calculated_sum_ratio.pdf')

#for n in ['name', 'code', 'E', 'M']:
#    df[n] = train_df[n]
#df = df.sort(columns='residual_squared')
#df.to_csv(os.path.join(dir,fname+'_residuals.csv'))



plt.show()
