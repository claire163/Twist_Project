"""
Some tools for working with chimera sequences
"""

import pandas as pd
from sys import exit

def contacting_terms (sample_space, contacts):
    """ Lists the possible contacts

    Parameters:
        sample_space (iterable): Each element in sample_space contains the possible
           amino acids at that position
        contacts (iterable): Each element in contacts pairs two positions that
           are considered to be in contact

    Returns:
        list: Each item in the list is a contact in the form ((pos1,aa1),(pos2,aa2))
    """
    contact_terms = []
    for contact in contacts:
        first_pos = contact[0]
        second_pos = contact[1]
        first_possibilities = set(sample_space[first_pos])
        second_possibilities = set(sample_space[second_pos])
        for aa1 in first_possibilities:
            for aa2 in second_possibilities:
                contact_terms.append(((first_pos,aa1),(second_pos,aa2)))
    return contact_terms

def make_contact_X (seqs, sample_space, contacts):
    contact_X = []
    contact_terms = contacting_terms(sample_space, contacts)
    for seq in seqs:
        cons = get_contacts(seq, contacts)
        inds = [contact_terms.index(c) for c in cons]
        X_row = [1 if i in inds else 0 for i in range(len(contact_terms))]
        contact_X.append(X_row)
    return contact_X, contact_terms

def make_sequence_X(seqs, sample_space):
    """ Make Xs for regression based on sequence elements."""
    sequence_X = []
    sequence_terms = [(i,t) for i,sp in enumerate(sample_space) for t in sp]
    for seq in seqs:
        this_terms = [(i,t) for i,t in enumerate(list(seq))]
        inds = [sequence_terms.index(s) for s in this_terms]
        X_row = [1 if i in inds else 0 for i in range(len(sequence_terms))]
        sequence_X.append(X_row)
    return sequence_X, sequence_terms

def make_X(seqs, sample_space, contacts):
    """ Make combined sequence/structure X. """
    seq_X, sequence_terms = make_sequence_X(seqs, sample_space)
    struct_X, contact_terms = make_contact_X(seqs, sample_space, contacts)
    X = [seq_X[i] + struct_X[i] for i in range(len(seqs))]
    terms = sequence_terms + contact_terms
    return X, terms

def get_contacts(seq, contacts):
    """
    Gets the contacts for seq.
    """
    cons = []
    for con in contacts:
        term = ((con[0],seq[con[0]]),(con[1],seq[con[1]]))
        cons.append(term)
    return cons

def load_assignments (assignments_file):
    assignments_line = [l for l in open(assignments_file).read().split('\n') if len(l)>0 and l[0]!='#']
    assignment = [ord(l.split('\t')[2]) - ord('A') for l in assignments_line if l.split('\t')[2] !='-']
    nodes_outputfile = [int(l.split('\t')[1])-1 for l in assignments_line if l.split('\t')[2] !='-'] # -1 because counting 0,1,2...
    return dict(zip(nodes_outputfile,assignment))

def make_sequence (code, assignments_dict, sample_space):
    ''' Returns the chimera amino acid sequence as a list'''
    seq = []
    for pos,aa in enumerate(sample_space):
        # Figure out which parent to use at that position
        if pos in assignments_dict:
            block = assignments_dict[pos] # the assigned block (based on pos)
            # the parent for that block in this particular chimera
            parent = int(code[block])
        else:
            parent = 0
        seq.append (aa[parent])
    return seq

def make_name_dict(dict_file):
    '''
    Makes the name dict from a spreadsheet
    '''
    name_df = pd.read_excel (dict_file)
    name_df['name'] = [s.lower() for s in name_df['name']]
    new_code = [zero_index(c) for c in name_df['code']]
    name_df['code'] = new_code
    name_dict = {a:b for a,b in zip(name_df['name'], new_code)}
    return name_dict

def zero_index (code):
    '''
    Takes a 1-indexed chimera code and zero-indexes it
    '''
    return ''.join ([str(int(x)-1) for x in str(code)])
