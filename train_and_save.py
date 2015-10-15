''' Reads in an ncr training set, builds a GPModel, and pickles the model

python2 train_and_save.py -k structure -t res.xlsx -n name

'''
import sys
sys.path.append ('/Users/seinchin/Documents/Caltech/Arnold Lab/Programming tools/GPModel')
import argparse, gpmodel, gpkernel, os, gptools
import cPickle as pickle
from chimera_tools import *
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import datetime

def main ():
    dir = os.path.dirname(__file__)
    parser = argparse.ArgumentParser()
    parser.add_argument('-k', '--kernel',required=True)
    parser.add_argument('-t', '--training',required=True)
    parser.add_argument('-n', '--name',required=True)
    parser.add_argument('-y', '--y_column',required=True)
    parser.add_argument('-p', '--plot', action='store_true')
    parser.add_argument('-d', '--drop',required=False, type=float)


    args = parser.parse_args()
    a_and_c = os.path.join(dir, 'alignment_and_contacts.pkl')
    sample_space, contacts = pickle.load(open(a_and_c))
    c_assignments_file = os.path.join(dir,'clibrary.output')
    n_assignments_file = os.path.join(dir,'nlibrary.output')
    dict_file = os.path.join(dir,'dict.xlsx')

    print 'Creating kernel...'
    if args.kernel == 'hamming':
        kern = gpkernel.HammingKernel()
    elif args.kernel == 'structure':
        kern = gpkernel.StructureKernel(contacts, sample_space)
    else:
        sys.exit ('Kernel type must be hamming or structure')


    print 'Creating training set...'
    # this reads in the chimeras for the training set and the file that relates
    # the abbreviations to which blocks are present
    # we need to shift all the parent numbers down by one to be consistent with
    # how they are recorded elsewhere
    with open (os.path.join(dir,args.training),'r') as f:
        train_df = pickle.load(f)


    if args.drop == None:
        # pull out the entries that don't have NaN in the y_column
        not_dropped = ~pd.isnull(train_df[args.y_column])
        Ys = train_df[not_dropped][args.y_column]
        Ys.index = train_df[not_dropped]['name']
        # make the X_seqs
        X_seqs = [list(seq) for seq in train_df[not_dropped]['sequence']]
        X_seqs = pd.DataFrame(X_seqs, index = Ys.index)

    else:
        Ys = train_df[args.y_column]
        Ys = Ys.fillna(args.drop)
        Ys.index = train_df['name']
        X_seqs = [list(seq) for seq in train_df['sequence']]
        X_seqs = pd.DataFrame( X_seqs, index = Ys.index)

    dt = datetime.date.today()
    name = str(dt) + '_' + args.name + '_' + args.y_column + '_' + args.kernel

    if args.plot:
        print 'Making LOO plot...'
        predicted, std = gptools.plot_LOO(X_seqs,Ys,kern, lab=args.y_column,
                        save_as=os.path.join(dir,name+'_LOO.pdf'))
        with open(os.path.join(dir,name+'_LOO.txt'),'w') as f:
            for i,n in enumerate(Ys.index):
                f.write (n+','+str(Ys[n])+','+str(predicted[i]))
                if std:
                    f.write(','+str(std[i]))
                f.write('\n')
        plt.show ()

    else:
        print 'Training model...'
        model = gpmodel.GPModel(X_seqs,Ys,kern,guesses=[100,10])
        print 'Pickling model...'
        with open(os.path.join(dir, name+ '_kernel.pkl'), 'wb') as f:
            pickle.dump(model, f)



if __name__ == "__main__":
    main()
