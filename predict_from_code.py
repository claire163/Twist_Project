import sys
sys.path.append ('/Users/seinchin/Documents/Caltech/Arnold Lab/Programming tools/GPModel')
import argparse, gpmodel, gpkernel, os, gptools
import cPickle as pickle
import chimera_tools
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import datetime
from itertools import product

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
        if i != len(seqs.index)-1:
            formatted += '\n'
    return formatted


def double_write(s, p=False, w=None):
    '''
    Print, write to file, or both
    '''
    if p:
        print s
    if w is not None:
        w.write(s+'\n')


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
with open(os.path.join(dir,args.model),'r') as m_file:
    model = pickle.load(m_file)


@profile
def main():
    if model.regr:
        head = 'code,mean,variance'
    else:
        head = 'code,pi'

    double_write(head, p=args.pr, w = out_file)

    # figure out which ones to predict
    predicted = False
    try:
        with open (os.path.join(dir,args.seq_file),'r') as f:
            # read and predict each line
            predicted = True
            pass
    except:
        pass

    if args.code != None:
        # read and predict the code
        double_write(formatted_predict(model, codes_to_seqs([args.code])),p=args.pr, w=out_file)
        predicted = True

    if not predicted:
        # generate and predict all possible chimeras
        for c in product([0,1,2],repeat=10):
            this = ''.join([str(i) for i in c])
            codes = ['n'+this, 'c'+this]
            cod = codes_to_seqs(codes)
            preds = formatted_predict(model,cod)
            double_write(preds,p=args.pr,w=out_file)
            #double_write(formatted_predict(model, codes_to_seqs(codes)),p=args.pr, w=out_file)


if __name__=="__main__":
    main()










