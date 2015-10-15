import pickle
import pandas as pd
import os
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import datetime
from sys import exit
from mpl_toolkits.mplot3d import Axes3D
import matplotlib

rc = {'lines.linewidth': 2,
      'axes.labelsize': 18,
      'axes.titlesize': 18,
      'axes.facecolor': 'DFDFE5'}
sns.set_context('notebook', rc=rc)
sns.set_style('darkgrid', rc=rc)

dt = datetime.date.today()

dir = os.path.dirname(__file__)
parser = argparse.ArgumentParser()
parser.add_argument('-n', '--name',required=True)
args = parser.parse_args()

with open (os.path.join(dir,args.name),'r') as f:
    df = pickle.load(f)

exp_inds = df['expression'] == 1
n_exp_inds = df['expression'] == -1



fig1 = plt.figure()
ax1 = sns.stripplot(x='exp_level', y='cschrimson_fraction', data=df, jitter=True, alpha=0.6)
plt.figure()
ax2 = sns.stripplot(x='exp_level', y='c1c2_fraction', data=df, jitter=True, alpha=0.6)
plt.figure()
ax3 = sns.stripplot(x='exp_level', y='cheriff_fraction', data=df, jitter=True, alpha=0.6)
plt.figure()
ax1 = sns.stripplot(x='loc_level', y='cschrimson_fraction', data=df[exp_inds], jitter=True, alpha=0.6)
plt.figure()
ax2 = sns.stripplot(x='loc_level', y='c1c2_fraction', data=df[exp_inds], jitter=True, alpha=0.6)
plt.figure()
ax3 = sns.stripplot(x='loc_level', y='cheriff_fraction', data=df[exp_inds], jitter=True, alpha=0.6)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for l,s in zip([-1,0,1,2],['red','green','violet','blue']):
    inds = df['exp_level'] == l
    #x = df[inds]['cschrimson_fraction']
    #y = df[inds]['c1c2_fraction']
    #z = df[inds]['cheriff_fraction']
    ax.scatter(df[inds]['c1c2_fraction'],df[inds]['cschrimson_fraction'],df[inds]['cheriff_fraction'],c=s)
figures=[manager.canvas.figure
         for manager in matplotlib._pylab_helpers.Gcf.get_all_fig_managers()]
for i, figure in enumerate(figures):
    figure.savefig(os.path.join(dir,'plots/parent_frac/'+str(dt)+'_figure%d.pdf' % i))
for y in ['mKate_mean', 'sum_ratio']:
    for ps in ['cschrimson_fraction', 'cheriff_fraction','c1c2_fraction']:
        figure = plt.figure()
        plt.plot(df[exp_inds][ps], df[exp_inds][y],'o')
        plt.xlabel(ps)
        plt.ylabel(y)
        figure.savefig(os.path.join(dir,'plots/parent_frac/'+str(dt)+'_'+y+'_'+ps+'.pdf'))




#plt.show()

# EM_bin_exp = plt.figure()
# EM_bin_exp_ax = EM_bin_exp.add_subplot(1,1,1)
# EM_bin_exp_ax.plot (df[exp_inds]['E'],df[exp_inds]['M'],'o')
# EM_bin_exp_ax.plot (df[n_exp_inds]['E'],df[n_exp_inds]['M'],'o')
# EM_bin_exp_ax.set_xlabel('E')
# EM_bin_exp_ax.set_ylabel('M')
# EM_bin_exp_ax.legend(['expressed','not expressed'],loc='best')
# x1,x2,y1,y2 = EM_bin_exp_ax.axis()
# EM_bin_exp_ax.set_xlim([-2,x2])
# EM_bin_exp_ax.set_ylim([-2,y2])
# EM_bin_exp_ax.set_title('Binary Expression with E and M')


# EM_exp = plt.figure()
# EM_exp_ax = EM_exp.add_subplot(1,1,1)
# for l in [-1,0,1,2]:
#     inds = df['exp_level'] == l
#     EM_exp_ax.plot(df[inds]['E'],df[inds]['M'],'o')

# EM_exp_ax.set_xlabel('E')
# EM_exp_ax.set_ylabel('M')
# EM_exp_ax.legend(['none','low','medium','high'],loc='best')
# x1,x2,y1,y2 = EM_exp_ax.axis()
# EM_exp_ax.set_xlim([-2,x2])
# EM_exp_ax.set_ylim([-2,y2])
# EM_exp_ax.set_title('Expression level with E and M')


# exp_clusters = plt.figure()
# ax = sns.stripplot(x='exp_level', y='mKate_mean', data=df[exp_inds], jitter=True, alpha=0.6)
# x1,x2,y1,y2 = ax.axis()
# ax.set_ylim([0,y2])
# ax.set_title('Expression Clusters')

# exp_E = plt.figure()
# exp_E_ax = sns.stripplot(x='exp_level', y='E', data=df, jitter=True, alpha=0.6)
# exp_E_ax.set_title('E for each expression level')

# exp_M = plt.figure()
# exp_M_ax = sns.stripplot(x='exp_level', y='M', data=df, jitter=True, alpha=0.6)
# exp_M_ax.set_title('M for each expression level')


# loc_clusters = plt.figure()
# ax2 = sns.stripplot(x='loc_level', y='sum_ratio', data=df[exp_inds], jitter=True, alpha=0.6)
# ax2.set_title ('Localization Clusters')

# loc_E = plt.figure()
# loc_E_ax = sns.stripplot(x='loc_level', y='E', data=df, jitter=True, alpha=0.6)
# loc_E_ax.set_title('E for each localization level')

# loc_M = plt.figure()
# loc_M_ax = sns.stripplot(x='loc_level', y='M', data=df, jitter=True, alpha=0.6)
# loc_M_ax.set_title('M for each localization level')

# EM_loc = plt.figure()
# EM_loc_ax = EM_loc.add_subplot(1,1,1)
# for l in [0,1,2]:
#     inds = df['loc_level'] == l
#     EM_loc_ax.plot(df[inds]['E'],df[inds]['M'],'o')

# EM_loc_ax.set_xlabel('E')
# EM_loc_ax.set_ylabel('M')
# EM_loc_ax.legend(['low','medium','high'],loc='best')
# EM_loc_ax.set_title('Localization level with E and M')

# EM_bin_exp.savefig(os.path.join(dir,'plots/'+str(dt)+'_EM_bin_exp.pdf'))
# EM_exp.savefig(os.path.join(dir,'plots/'+str(dt)+'_EM_exp.pdf'))
# exp_clusters.savefig(os.path.join(dir,'plots/'+str(dt)+'_exp_clusters.pdf'))
# exp_E.savefig(os.path.join(dir,'plots/'+str(dt)+'_exp_E.pdf'))
# exp_M.savefig(os.path.join(dir,'plots/'+str(dt)+'_exp_M.pdf'))
# loc_clusters.savefig(os.path.join(dir,'plots/'+str(dt)+'_loc_clusters.pdf'))
# loc_E.savefig(os.path.join(dir,'plots/'+str(dt)+'_loc_E.pdf'))
# loc_M.savefig(os.path.join(dir,'plots/'+str(dt)+'_loc_M.pdf'))
# EM_loc.savefig(os.path.join(dir,'plots/'+str(dt)+'_EM_loc.pdf'))



