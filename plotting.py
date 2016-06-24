import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import argparse
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
parser.add_argument('-he', '--header', type=int, default=0)
parser.add_argument('-r', '--reals', required=False)
parser.add_argument('-he2', '--header_2', type=int, required=False)
parser.add_argument('-l', '--limit', type=float, required=False)



args = parser.parse_args()

parents = {'cschrimson':'red', 'c1c2':'green', 'cheriff':'blue'}
parent_names =  ['cschrimson', 'c1c2', 'cheriff']
df = pd.read_csv(args.data, skiprows=args.header, comment='#')
# regression results
if len(df.columns) == 4:
    plt.plot(df['y'], df['mu'], 'o', ms=6, color='grey', alpha=0.7)
    if (any(df['name'] == 'cschrimson')):
        for p in parent_names:
            plt.plot(df[df['name']==p]['y'], df[df['name']==p]['mu'],
                     'o', ms=8,color=parents[p], alpha=0.8)
        lg = plt.legend(('chimeras', 'CsChrimR', 'C1C2', 'CheRiff'),
                        loc='best')
    plt.xlabel('measured ' + args.name)
    plt.ylabel('predicted ' + args.name)
#     plt.gca().set_xlim([-11.0,0.5])
#     plt.gca().set_ylim([plt.gca().get_ylim()[0], 0.5])


elif len(df.columns) == 3:
    df_2 = pd.read_csv(args.reals, skiprows=args.header_2, comment='#')
    df['real'] = [df_2[df_2['name']==n]['y'] for n in df['name']]
    plt.plot(df['real'], df['pi'], 'o', ms=6, color='grey', alpha=0.7)
    if (any(df['name'] == 'cschrimson')):
        for p in parent_names:
            plt.plot(df[df['name']==p]['real'], df[df['name']==p]['pi'],
                     'o', ms=8, color=parents[p], alpha=0.8)
        lg = plt.legend(('chimeras', 'CsChrimR', 'C1C2', 'CheRiff'),
                        loc='best')
        lowest_parent = min(df_2[df_2['name'].isin(parent_names)]['y'])
        xlims = plt.gca().get_xlim()
        ylims = plt.gca().set_ylim([0, 1])
        plt.axvspan(xlims[0], lowest_parent, facecolor='grey', alpha=0.2)
    plt.xlabel('measured ' + args.name)
    plt.ylabel('predicted probability of ' + args.name)

elif len(df.columns) == 6:
    plt.plot(df['real'], df['pi'], 'o', color='grey', alpha=0.7)
    xlims = plt.gca().set_xlim([-3,1])
    ylims = plt.gca().set_ylim([0, 1])
    plt.axvspan(xlims[0], args.limit, facecolor='grey', alpha=0.2)
    plt.xlabel('measured ' + args.name)
    plt.ylabel('predicted probability of ' + args.name)

if args.write is not None:
    plt.savefig(args.write, bbox_inches='tight')
plt.show()




