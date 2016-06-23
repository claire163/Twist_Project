import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
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



args = parser.parse_args()

parents = {'cschrimson':'red', 'c1c2':'green', 'cheriff':'blue'}
df = pd.read_csv(args.data, skiprows=args.header, comment='#')
# regression results
if len(df.columns) == 4:
    plt.plot(df['y'], df['mu'], 'o', ms=6, color='grey', alpha=0.7)
    if (any(df['name'] == 'cschrimson')):
        for p in ['cschrimson', 'c1c2', 'cheriff']:
            plt.plot(df[df['name']==p]['y'], df[df['name']==p]['mu'],
                     'o', ms=8,color=parents[p], alpha=0.8)
        lg = plt.legend(('chimeras', 'CsChrimR', 'C1C2', 'CheRiff'),
                        loc='best')

    plt.margins(0.02)
    plt.xlabel('measured ' + args.name)
    plt.ylabel('predicted ' + args.name)
    if args.write is not None:
        plt.savefig(args.write, bbox_inches='tight')
    plt.show()


