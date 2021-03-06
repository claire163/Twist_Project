'''
Reads in an ncr training set, builds a GPModel, and pickles the model

python2 train_and_save.py -k structure -t res.xlsx -n name

'''
import sys
sys.path.append('/Users/kevinyang/Documents/Projects/GPModel')
sys.path.append ('/Users/seinchin/Documents/Caltech/Arnold Lab/Programming tools/GPModel')
import argparse, gpmodel, gpkernel, os, gptools, gpmean
import dill as pickle
from chimera_tools import *
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import datetime
from scipy import stats
from collections import namedtuple
from sklearn import metrics, linear_model


def main ():
    dir = os.path.dirname(__file__)
    parser = argparse.ArgumentParser()
    parser.add_argument('-k', '--kernel',required=True, nargs='*')
    parser.add_argument('-a', '--alpha', action='store_true')
    parser.add_argument('-t', '--training',required=True)
    parser.add_argument('-n', '--name',default=None)
    parser.add_argument('-y', '--y_column',required=True)
    parser.add_argument('-s', '--std_column', required=False)
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
    kerns = []
    for k in args.kernel:
        if k in ['hamming', 'n_code', 'c_code', 'subblock']:
            kern = gpkernel.HammingKernel()
        elif k == 'WHamming':
            kern = gpkernel.WeightedHammingKernel()
        elif k == 'structure':
            kern = gpkernel.StructureKernel(contacts)
        elif k == '32Structure':
            kern = gpkernel.StructureMaternKernel(contacts, '3/2')
        elif k == '52Structure':
            kern = gpkernel.StructureMaternKernel(contacts, '5/2')
        elif k == '32Hamming':
            kern = gpkernel.HammingMaternKernel('3/2')
        elif k == '52Hamming':
            kern = gpkernel.HammingMaternKernel('5/2')
        elif k == 'SEStructure':
            kern = gpkernel.StructureSEKernel(contacts)
        elif k in ['SEHamming', 'SEn_code', 'SEc_code', 'SEsubblock']:
            kern = gpkernel.HammingSEKernel()
        else:
            raise ValueError('Invalid kernel type: ' + k)
        kerns.append(kern)
    if len(kerns) == 1:
        kern = kerns[0]
    else:
        kern = gpkernel.SumKernel(kerns)


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

    # make the variances
    if args.std_column is not None:
        variances = train_df[not_dropped][args.std_column] ** 2
        variances.index = Ys.index
    else:
        variances = None

    # make the X_seqs
    if args.kernel in ['SEn_code', 'SEc_code', 'n_code', 'c_code']:
        s = 'code'
    elif args.kernel in ['subblock', 'SEsubblock']:
        s = 'subblock'
    else:
        s = 'sequence'


    X_seqs = [list(seq) for seq in train_df[not_dropped][s]]
    X_seqs = pd.DataFrame(X_seqs, index=Ys.index)


    print 'Training model...'
    if args.alpha:
        clf = linear_model.BayesianRidge
        mf = gpmean.StructureSequenceMean(sample_space, contacts,
                                          clf)
        model = gpmodel.GPModel(kern,
                                guesses=args.guess,
                                mean_func=mf,
                                objective=args.objective)
    else:
        model = gpmodel.GPModel(kern,
                                guesses=args.guess,
                                objective=args.objective)
    model.fit(X_seqs, Ys, variances=variances)
    print model.hypers
    print 'log_ML = %f' %-model.ML
    try:
        print 'log_LOO_P = %f' %-model.log_p
    except:
        pass

    if model.regr:
        LOOs = model.LOO_res (model.hypers, add_mean=True)
        actual = Ys
        predicted = LOOs['mu']
        var = LOOs['v']
        r1 = stats.rankdata(actual)
        r2 = stats.rankdata(predicted)
        tau = stats.kendalltau(r1, r2).correlation
        R = np.corrcoef(actual, predicted)[0,1]
        print 'tau = %.4f' %tau
        print 'R = %.4f' %R
        try:
            print 'avg var = %0.4f' %np.mean(model.variances)
        except TypeError:
            pass
        if args.alpha:
            print 'R_lin = %.4f' %np.corrcoef(actual,
                                              model.mean_func.means)[0,1]
    else:
        preds = model.predicts(model.X_seqs)
        preds = [p[0] for p in preds]
        fpr, tpr, _ = metrics.roc_curve(model.Y, preds)
        auc = metrics.auc(fpr,tpr)
        print 'AUC = %.4f' %auc

    if args.name is not None:
        print 'Pickling model...'
        dt = datetime.date.today()
        model_file = args.training.split('/')[0] + '/models/'
        name = model_file + '_'.join([str(dt), args.name, args.y_column,
                                      '_'.join(args.kernel)])
        model.dump(os.path.join(dir, name + '.pkl'))
        if model.regr:
            with open(os.path.join(dir, name + '_LOO.txt'),'w') as f:
                f.write('name,mu,var,y\n')
                for i,n in enumerate(Ys.index):
                    f.write (n+','+str(predicted[n]))
                    f.write(','+str(var[n])+','+str(Ys.loc[n]))
                    f.write('\n')
        else:
            predicted, actual, _ = gptools.cv(model.X_seqs, model.Y, model,
                                              len(X_seqs) - 1)
            fpr, tpr, _ = metrics.roc_curve(model.Y, predicted)
            auc = metrics.auc(fpr,tpr)
            print 'LOO AUC = %.4f' %auc
            with open (name + '_LOO.txt', 'w') as f:
                f.write ('name,y,pi\n')
                for n, r, p in zip (model.Y.index, model.Y, predicted):
                    f.write (n + ',' + str(r) + ',' + str(p) + '\n')
    if args.plot:
        print 'Making LOO plot...'
        if model.regr:
            gptools.plot_predictions(actual, predicted,
                                     label=args.y_column)
            parents = ['c1c2', 'cschrimson', 'cheriff']
            colors = ['green', 'red', 'blue']
            for p, c in zip (parents, colors):
                a = actual[p]
                pr = predicted[p]
                plt.plot (a, pr, '.', color=c)
            plt.margins(0.02)
            if args.name is not None:
                plt.savefig(name + '_LOO.pdf')

        else:
            auc = gptools.plot_ROC(model.Y, preds)
            if args.name is not None:
                plt.savefig(name+'_ROC.pdf')
        if args.name is None:
            plt.show()
    print 'Saving model results...'
    k = ' '.join(args.kernel)
    if model.regr:
        save_me = [k, args.y_column, str(-model.ML),
                   str(-model.log_p), str(R), str(tau),
                   '', ' '.join(str(model.hypers).split(',')), 'no',
                   str(args.name is not None)]
    else:
        save_me = [k, args.y_column, str(-model.ML),
                   '', '', '',str(auc),
                   ' '.join(str(model.hypers).split(',')), 'no',
                   str(args.name is not None)]
    with open(args.training.split('/')[0] + '/models.csv', 'r') as f:
        lines = [','.join(ell.split(',')[0:-2]) for ell in f.readlines()]
    if not ','.join(save_me[0:-2]) in lines:
        with open(args.training.split('/')[0] + '/models.csv', 'a') as f:
            f.write('\n' + ','.join(save_me))

if __name__ == "__main__":
    main()
