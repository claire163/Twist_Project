""" read in LOO data and plot it"""
import matplotlib.pyplot as plt
import pandas as pd

import sys
sys.path.append ('/Users/seinchin/Documents/Caltech/Arnold Lab/Programming tools/GPModel')
import gpmodel, gpkernel, os, gptools, pickle
import numpy as np
import datetime

import seaborn as sns
sns.set_context('notebook', font_scale=1.5, rc={'lines.linewidth': 2.5})
sns.set_style('darkgrid', {'axes.facecolor': '(0.875, 0.875, 0.9)'})

dir = os.path.dirname(__file__)

with open (os.path.join(dir,'2015-10-20_all_res.pkl'),'r') as f:
    df = pickle.load(f)

exp_inds = df['expression'] == 1
n_exp_inds = df['expression'] == -1


c_df = pd.read_csv('c_chimeras_data.txt',header=None,delimiter=' ')
n_df=pd.read_csv('n_chimeras_data.txt',delimiter=' ',header=None)
c_df.columns=['code','E','M','Sequence']
n_df.columns = c_df.columns

plt.figure(figsize=(15,10))
plt.plot(c_df['E'],c_df['M'],'.',alpha=.9)
plt.plot(n_df['E'],n_df['M'],'.',alpha=.4)

inds = df['exp_level'] < 1
plt.plot(df[inds]['E'],df[inds]['M'],'o')
inds = df['exp_level'] >= 1
plt.plot(df[inds]['E'],df[inds]['M'],'o')


# x1,x2,y1,y2 = EM_exp_ax.axis()
# EM_exp_ax.set_xlim([-2,x2])
# EM_exp_ax.set_ylim([-2,y2])
plt.title('Expression level with E and M')


plt.xlabel('E')
plt.ylabel('M')
plt.legend(['contiguous library','non-contiguous library',
           'poor expression','good expression'],loc='best')
plt.savefig('2015-10-22_E_M_expression2.pdf')
#plt.show()

#fname = 'LOO_results/2015-10-14__expression_hamming_LOO.txt'
#with open (os.path.join(dir,fname),'r') as f:
#    c_df = pd.read_csv(f)

# df.columns = ['name','actual','predicted']

# gptools.plot_ROC (df['actual'], df['predicted'], file_name='2015-10-14_exp_hamming_kernel_LOO_ROC.pdf',title='Expression Hamming LOO')
# plt.show()



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
