import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import argparse
import pickle
from sys import exit
sns.set_style('whitegrid')
sns.set_context('paper')
# Plot adjustments:
plt.rcParams.update({'ytick.labelsize': 14})
plt.rcParams.update({'xtick.labelsize': 14})
plt.rcParams.update({'axes.labelsize': 16})
plt.rcParams.update({'legend.fontsize': 14})


parser = argparse.ArgumentParser()
parser.add_argument('-d', '--data', required=True)
parser.add_argument('-n', '--name', required=True)
parser.add_argument('-w', '--write', required=False)
args = parser.parse_args()


df = pd.read_csv(args.data)
parent_names = ['n0000000000', 'n1111111111', 'n2222222222']
parent_colors = ['blue', 'green', 'red']
parents_dict = {n:c for n,c in zip(parent_names, parent_colors)}
with open('2016-06-22/props.pkl', 'rb') as f:
    trained_df = pickle.load(f, encoding='latin1')
trained_codes = [n[0] + str(c) for n,c in zip(trained_df['name'], trained_df['code'])]


# regression results
if len(df.columns) == 3:
    df = df.sort_values('mu')
    df['std'] = np.sqrt(df['variance'])
    df['rank'] = range(len(df))
    plt.plot(df['rank'], df['mu'], color='grey', alpha=0.7)
    trained = df[df['code'].isin(trained_codes)]
#     plt.errorbar(trained['rank'], trained['mu'],
#                  yerr=[s for s in trained['std']],
#                 fmt='o', color='grey', alpha=0.8)
    for p in parent_names[::-1]:
        plt.errorbar(df[df['code']==p]['rank'], df[df['code']==p]['mu'],
                     yerr=df[df['code']==p]['std'], fmt='o',
                     color=parents_dict[p], alpha=0.8, ms=6)
    plt.fill_between(df['rank'], df['mu'] + df['std'], df['mu'] - df['std'],
                     facecolor='grey', alpha=0.3)
    plt.gca().set_xlim([plt.gca().get_xlim()[0] - 1000,
                        plt.gca().get_xlim()[1]])
    plt.gca().get_xaxis().set_ticks([])
    plt.ylabel('predicted ' + args.name)
    lg = plt.legend(('chimeras', 'CsChrimR', 'C1C2', 'CheRiff'),
                    loc='best')
#     plt.gca().set_xlim([-11.0,0.5])
#     plt.gca().set_ylim([plt.gca().get_ylim()[0], 0.5])




if args.write is not None:
    plt.savefig(args.write, bbox_inches='tight')
plt.show()




