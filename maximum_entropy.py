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
    parser.add_argument('-n', '--n', required=True, type=int)
    parser.add_argument('-p', '--pr', action='store_true')
    parser.add_argument('-w', '--out_file', required=False)
    args = parser.parse_args()



    # load the sequences
    print 'Attempting to load all possible sequences...'
    try:
        sequences = pd.read_csv('all_chimeras.txt')
        sequences.index = xrange(len(sequences))
        codes = sequences.iloc[:,0].values
        sequences = sequences.drop(sequences.columns[[0]], axis=1)
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
    kernel = model.kern
    vn = model.hypers[0]
    hypers = model.hypers[1::]
    observations = model.X_seqs
    # remove observations from sequences
    obs_list = [''.join(x) for _, x in observations.iterrows()]
    seq_list = [''.join(x) for _, x in sequences.iterrows()]
    seq_list = pd.DataFrame(seq_list, columns=['seq'])
    keep = ~seq_list['seq'].isin(obs_list)
    sequences = sequences[keep]
    print 'Loading probabilities...'
    predictions = pd.read_csv(args.prediction_file)
    probabilities = predictions[keep][['pi']].values
    print 'Maximizing entropy...'
    ent = gpentropy.GPEntropy(model=model)
    selected, H, inds = ent.maximize_expected_entropy(sequences, probabilities, args.n)
    print selected
    print H
    print codes[inds]
    if args.out_file is not None:
        with open(args.outfile, 'w') as f:
            f.write('model:')
            f.write(args.model + '\n')
            f.write('probabilities:' + args.prediction_file + '\n')
            f.write('H=%f' %H)
            for se, i in zip(selected, inds):
                f.write('\n')
                f.write(codes[i])
                f.write(',')
                f.write(se)


if __name__=="__main__":
    main()