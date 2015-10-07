"""
Some tools for working with chimera sequences
"""

import pandas as pd

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
