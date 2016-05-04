import sys, os
import dill as pickle
import pandas as pd
import numpy as np
sys.path.append('/Users/kevinyang/Documents/Projects/GPModel')
sys.path.append ('/Users/seinchin/Documents/Caltech/Arnold Lab/Programming tools/GPModel')
import argparse, gpmodel, gpkernel, gpentropy, chimera_tools, itertools

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
    return chimera_tools.make_sequence(code[1:], my_dict[code[0]], sample_space)

dir = os.path.dirname(__file__)
c_assignments_file = os.path.join(dir,'clibrary.output')
n_assignments_file = os.path.join(dir,'nlibrary.output')
c_assignments_dict = chimera_tools.load_assignments(c_assignments_file)
n_assignments_dict = chimera_tools.load_assignments(n_assignments_file)
a_and_c = os.path.join(dir, 'alignment_and_contacts.pkl')
sample_space, contacts = pickle.load(open(a_and_c))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--model', required=True)
    parser.add_argument('-l', '--prediction_file', required=True)
    parser.add_argument('-p', '--pr', action='store_true')
    parser.add_argument('-w', '--out_file', required=True)
    args = parser.parse_args()



    # load the sequences
    print 'Attempting to load all possible sequences...'
    try:
        sequences = pd.read_csv('all_chimeras.txt')
        sequences.index = xrange(len(sequences))
        print '\tSuccess!'
    except:
        print 'Generating all possible sequences'
        # generate all possible chimeras
        # load the assignment dicts for contiguous and non-contiguous
        sequences = pd.DataFrame()
        for c in itertools.product([0,1,2],repeat=10):
            this = ''.join([str(i) for i in c])
            codes = ['n'+this, 'c'+this]
            sequences = pd.concat([sequences, codes_to_seqs(codes)])
        sequences.to_csv('all_chimeras.txt')

    #load the model
    print 'Loading the model...'
    model = gpmodel.GPModel.load(args.model)


    # load the predictions
    predictions = pd.read_csv(args.prediction_file)

if __name__=="__main__":
    main()