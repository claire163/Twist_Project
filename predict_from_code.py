import sys
sys.path.append('/Users/kevinyang/Documents/Projects/GPModel')
sys.path.append ('/Users/seinchin/Documents/Caltech/Arnold Lab/Programming tools/GPModel')
import argparse, gpmodel, gpkernel, os, gptools
import cPickle as pickle
import chimera_tools
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import datetime
import itertools
from sys import exit

""" Here, 'code' refers to 'n' or 'c' + (10 digits)"""

def codes_to_seqs(codes):
    """
    Take a list of chimera codes and return a DataFrame of sequences
    """
    if type(codes) is str:
        sys.exit('codes must be a list, not a string')
    return pd.DataFrame([code_to_seq(c) for c in codes], index=codes)

def code_to_seq(code):
    """
    Take a chimera code and return the sequence as a list
    """
    my_dict = {'n':n_assignments_dict,'c':c_assignments_dict}
    return chimera_tools.make_sequence(code[1:],my_dict[code[0]],sample_space)

def formatted_predict(model, seqs):
    '''
    Make predictions and return a string where each row is
    ind,mean_or_pi,var
    '''
    preds = model.predicts(seqs)
    formatted = ''
    for i, ind in enumerate(seqs.index):
        formatted += ','.join([ind] + [str(p) for p in preds[i]])
        formatted += '\n'
    return formatted

def double_write(s, p=False, w=None):
    '''
    Print, write to file, or both
    '''
    if p:
        print s
    if w is not None:
        w.write(s)

def get_all_sequences():
    try:
        print 'Attempting to load all possible sequences...'
        sequences = pd.read_csv('all_chimeras.txt')
        sequences.index = sequences['Unnamed: 0']
        sequences = sequences.drop(sequences.columns[[0]], axis=1)
        print '\tSuccess!'
    except:
        print 'Generating all possible sequences...'
        # generate all possible chimeras
        # load the assignment dicts for contiguous and non-contiguous
        sequences = pd.DataFrame()
        for c in itertools.product([0,1,2],repeat=10):
            this = ''.join([str(i) for i in c])
            codes = ['n'+this, 'c'+this]
            sequences = pd.concat([sequences, codes_to_seqs(codes)])
        sequences.to_csv('all_chimeras.txt')
    return sequences

dir = os.path.dirname(__file__)
parser = argparse.ArgumentParser()
parser.add_argument('-m', '--model',required=True)
parser.add_argument('-s', '--seq_file',required=False)
parser.add_argument('-c', '--code',required=False)
parser.add_argument('-p', '--pr', action='store_true')
parser.add_argument('-w', '--write', required=False)

args = parser.parse_args()

if args.write:
    out_file = open(os.path.join(dir, args.write), 'w')
else:
    out_file = None

# load the assignment dicts for contiguous and non-contiguous
c_assignments_file = os.path.join(dir,'clibrary.output')
n_assignments_file = os.path.join(dir,'nlibrary.output')
c_assignments_dict = chimera_tools.load_assignments(c_assignments_file)
n_assignments_dict = chimera_tools.load_assignments(n_assignments_file)
a_and_c = os.path.join(dir, 'alignment_and_contacts.pkl')
sample_space, contacts = pickle.load(open(a_and_c))

#load the model
print 'Loading the model...'
model = gpmodel.GPModel.load(args.model)
try:
    with open(args.model.split('.pkl')[0] + '_terms.pkl', 'rb') as f:
        terms = pickle.load(f)
        has_terms = True
except IOError:
    has_terms = False

def main():
    if model.regr:
        head = 'code,mu,variance\n'
    else:
        head = 'code,pi,f_bar,variance\n'

    if has_terms:
        global terms

    double_write(head, p=args.pr, w = out_file)

    # figure out which ones to predict
    predicted = False
    try:
        with open (os.path.join(dir,args.seq_file),'r') as f:
            seqs = pd.read_csv(f, header=None)
        seqs.columns = ['name', 'sequence']
        if has_terms:
            X, terms = make_X(seqs['sequence'], sample_space, contacts,
                                 terms=terms, collapse=False)
            seqs = pd.DataFrame(X, index=seqs.name)
        else:
            seqs = pd.DataFrame([[s for s in se] for se in seqs.sequence],
                            index=seqs.name)
        double_write(formatted_predict(model, seqs), p=args.pr, w=out_file)
        predicted = True
    except:
        pass
    if args.code != None:
        # read and predict the code
        double_write(formatted_predict(model, codes_to_seqs([args.code])),
                     p=args.pr, w=out_file)
        predicted = True

    if not predicted:
        if has_terms:
            print 'Generating all possible Xs...'
            sequences = get_all_sequences()
            seqs = [''.join(s) for _, s in sequences.iterrows()]
            seqs = pd.DataFrame(seqs, index=sequences.index,
                                    columns=['sequence'])
            n_splits = 1000
            n_per = len(seqs.index) / n_splits
            inds = [seqs.index[n * n_per: (n+1) * n_per]
                    for n in range(n_splits)]
            inds.append(seqs.index[n_splits * n_per::])
            all_Xs = []
            for i,ind in enumerate(inds):
                X, terms = chimera_tools.make_X(seqs.loc[ind]['sequence'],
                                                sample_space,
                                                contacts, terms=terms,
                                                collapse=False)
                all_Xs.append(pd.DataFrame(X, index=ind))
            sequences = pd.concat(all_Xs)
        else:
            sequences = get_all_sequences()

        for i in sequences.index:
            seq = sequences.loc[[i]]
            preds = formatted_predict(model,seq)
            double_write(preds,p=args.pr,w=out_file)
    if args.write:
        out_file.close()

if __name__=="__main__":
    main()
