''' Reads in an ncr training set, builds a GPModel, and pickles the model

python2 train_and_save.py -k structure -t res.xlsx -n name

'''
import sys
sys.path.append ('/Users/seinchin/Documents/Caltech/Arnold Lab/Programming tools/GPModel')
import argparse, gpmodel, gpkernel, os, pickle, gptools
from chimera_tools import *
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def main ():
    dir = os.path.dirname(__file__)
    parser = argparse.ArgumentParser()
    parser.add_argument('-k', '--kernel',required=True)
    parser.add_argument('-t', '--training',required=True)
    parser.add_argument('-n', '--name',required=True)
    parser.add_argument('-y', '--y_column',required=True)

    args = parser.parse_args()
    a_and_c = os.path.join(dir, 'alignment_and_contacts.pkl')
    sample_space, contacts = pickle.load(open(a_and_c))
    training_chimeras = os.path.join(dir, args.training)
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
    train_df = pd.read_excel(training_chimeras)
    train_df = pd.DataFrame.dropna(train_df) # drop rows that have not been tested
    train_df['name'] = [s.lower() for s in train_df['name']]
    Ys = pd.Series([t for t in train_df[args.y_column]], index=train_df['name'])
    # make binary data 1/-1 instead of 1/0
    if sum([1 if y==1 or y==0 else 0 for y in Ys]) == len(Ys):
        Ys = Ys*2-1

    # make the name_dict from dict.xlsx
    name_dict = make_name_dict(dict_file)


    ## load the assignments files
    # put this info in a dictionary because it's easier
    c_assignments_dict = load_assignments(c_assignments_file)
    n_assignments_dict = load_assignments(n_assignments_file)

    # for each member of the training set, we want to generate a sequence
    X_seqs = []
    for name in train_df['name']:
        if name[0] in set(['n', 'N']):
            X_seqs.append(make_sequence(name_dict[name], n_assignments_dict, sample_space))
        else:
            X_seqs.append(make_sequence(name_dict[name], c_assignments_dict, sample_space))
    X_seqs = pd.DataFrame(X_seqs, index = train_df['name'])
    gptools.plot_LOO(X_seqs,Ys,kern)#,save_as=os.path.join(dir,'ham_class_test.png'))
    plt.show ()
    exit('')

    print 'Training model...'
    model = gpmodel.GPModel(X_seqs,Ys,kern,guesses=[100,10])
    print 'Pickling model...'
    with open(os.path.join(dir, args.name + args.kernel + '_kernel.pkl'), 'wb') as f:
        pickle.dump(model, f)



if __name__ == "__main__":
    main()
