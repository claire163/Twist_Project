'''
Reads in an ncr training set, builds a GPModel, and pickles the model

python2 train_and_save.py -k structure -t res.xlsx -n name

'''
import sys
sys.path.append ('/Users/seinchin/Documents/Caltech/Arnold Lab/Programming tools/GPModel')
import argparse, gpmodel, gpkernel, os, gptools
import dill as pickle
from chimera_tools import *
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import datetime
from scipy import stats
from collections import namedtuple

def main ():
    dir = os.path.dirname(__file__)
    parser = argparse.ArgumentParser()
    parser.add_argument('-k', '--kernel',required=True)
    parser.add_argument('-t', '--training',required=True)
    parser.add_argument('-n', '--name',default=None)
    parser.add_argument('-y', '--y_column',required=True)
    parser.add_argument('-p', '--plot', action='store_true')
    parser.add_argument('-d', '--drop',required=False, type=float)
    parser.add_argument('-g', '--guess', nargs='*', required=False)
    parser.add_argument('-o', '--objective', required=False, default='log_ML')


    args = parser.parse_args()
    a_and_c = os.path.join(dir, 'alignment_and_contacts.pkl')
    sample_space, contacts = pickle.load(open(a_and_c))
    c_assignments_file = os.path.join(dir,'clibrary.output')
    n_assignments_file = os.path.join(dir,'nlibrary.output')
    dict_file = os.path.join(dir,'dict.xlsx')

    print 'Creating kernel...'
    if args.kernel in ['hamming', 'n_code', 'c_code', 'subblock']:
        kern = gpkernel.HammingKernel()
    elif args.kernel == 'WHamming':
        kern = gpkernel.WeightedHammingKernel()
    elif args.kernel == 'structure':
        kern = gpkernel.StructureKernel(contacts)
    elif args.kernel == 'SEStructure':
        kern = gpkernel.StructureSEKernel(contacts)
    elif args.kernel in ['SEHamming', 'SEn_code', 'SEc_code', 'SEsubblock']:
        kern = gpkernel.HammingSEKernel()
    else:
        sys.exit ('Invalid kernel type')


    print 'Creating training set...'
    # this reads in the chimeras for the training set and the file that relates
    # the abbreviations to which blocks are present
    # we need to shift all the parent numbers down by one to be consistent with
    # how they are recorded elsewhere
    with open (os.path.join(dir,args.training),'r') as f:
        train_df = pickle.load(f)


    # Make the Ys
    if args.drop == None:
        # pull out the entries that don't have NaN in the y_column
        not_dropped = ~pd.isnull(train_df[args.y_column])
    else:
        not_dropped = [True for _ in train_df.index]
    if args.kernel == 'n_code' or args.kernel == 'SEn_code':
        in_lib = [n[0] == 'n' or n in ['c1c2','cheriff','cschrimson']\
                 for n in train_df['name']]
        not_dropped = [d and i for d,i in zip(not_dropped, in_lib)]
    elif args.kernel == 'c_code' or args.kernel == 'SEc_code':
        in_lib = [n[0] == 'c' or n in ['c1c2','cheriff','cschrimson']\
                 for n in train_df['name']]
        not_dropped = [d and i for d,i in zip(not_dropped, in_lib)]
    not_dropped = pd.Series(not_dropped, index=train_df.index)
    Ys = train_df[not_dropped][args.y_column]
    if args.drop is not None:
        Ys = Ys.fillna(args.drop)
    Ys.index = train_df[not_dropped]['name']

    # make the X_seqs
    if args.kernel in ['SEn_code', 'SEc_code', 'n_code', 'c_code']:
        s = 'code'
    elif args.kernel in ['subblock', 'SEsubblock']:
        s = 'subblock'
    else:
        s = 'sequence'


    X_seqs = [list(seq) for seq in train_df[not_dropped][s]]
    X_seqs = pd.DataFrame(X_seqs, index = Ys.index)


    print 'Training model...'
    model = gpmodel.GPModel(kern,
                            guesses=args.guess,
                            objective=args.objective)
    model.fit(X_seqs, Ys)
    print model.hypers
    print '-log_ML = %f' %model.ML
    try:
        print '-log_LOO_P = %f' %model.log_p
    except:
        pass

    if args.name is not None:
        print 'Pickling model...'
        dt = datetime.date.today()
        name = str(dt) + '_' + args.name + '_' + args.y_column + '_' + args.kernel
        model.dump(os.path.join(dir, name + '_dict.pkl'))


    if args.plot:
        print 'Making LOO plot...'

        if model.regr:
            LOOs = model.LOO_res (model.hypers)
            predicted = model.unnormalize(LOOs['mu'])
            var = model.std**2*LOOs['v']
            actual = model.Y
            r1 = stats.rankdata(actual)
            r2 = stats.rankdata(predicted)
            print 'tau = %.4f' %stats.kendalltau(r1, r2).correlation

            gptools.plot_predictions(actual, predicted,
                                     label=args.y_column)
            parents = ['c1c2', 'cschrimson', 'cheriff']
            colors = ['green', 'red', 'blue']
            for p, c in zip (parents, colors):
                a = actual[p]
                pr = predicted[p]
                plt.plot (a, pr, '.', color=c)
            plt.title(args.kernel + ' model')
            print 'R = %.4f' %np.corrcoef(model.normed_Y, LOOs['mu'])[0,1]
            if args.name is not None:
#                 plt.savefig('plots/' + args.kernel + \
#                             '_model_' + args.y_column + '_LOO.pdf')
                plt.savefig('plots/'+name+'_LOO.pdf')
                with open(os.path.join(dir,'LOO_results/'+name+'_LOO.txt'),'w') as f:
                    f.write('name,mu,var,y\n')
                    for i,n in enumerate(Ys.index):
                        f.write (n+','+str(predicted[n]))
                        f.write(','+str(var[n])+','+str(Ys.loc[n]))
                        f.write('\n')

        else:
            preds = model.predicts(model.X_seqs)
            auc = gptools.plot_ROC(model.Y, preds)
            print 'AUC = %.4f' %auc

            if args.name is not None:
                with open ('LOO_results/' + name + '_res.txt', 'w') as f:
                    f.write ('name, expression, pi\n')
                    for n, r, p in zip (model.Y.index, model.Y, preds):
                        f.write (n + ',' + str(r) + ',' + str(p) + '\n')
                plt.title ('ROC for ' + args.y_column + \
                           ' with ' + args.kernel + ' model')
                plt.savefig('plots/'+name+'_ROC.pdf')

        plt.show()


if __name__ == "__main__":
    main()
