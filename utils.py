
from IPython.display import IFrame
import numpy as np

def calc_MCC(prediction, ref):
    TP = 0
    FP = 0
    FN = 0
    TN = 0

    for i in range(len(prediction)):
        if ref[i] != -1:                # reference is paired
            if prediction[i] == ref[i]: # correct pair predicted
                TP += 1
            elif prediction[i] == -1:   # predicted as unpaired
                FN += 1
            else:                       # incorrect pair predicted
                FP += 1 
        elif ref[i] == -1:              # reference is unpaired
            if prediction[i] == ref[i]: # predicted as unpaired
                TN += 1
            else:                       # predicted is paired to something
                FP += 1

    MCC = ((TP * TN) - (FP * FN)) / np.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    if np.isnan(MCC):
        return(-1)

    return MCC

# Turn a db string into a dict
def dict_dot_bracket(db):
    open_stack = []
    pairs = {}
    for i, c in enumerate(db):
        if c == '(':
            open_stack.append(i)
        elif c == ')':
            j = open_stack.pop()
            pairs[i+1] = j+1 # ViennaRNA appears to be 1-indexed!
            pairs[j+1] = i+1

    return(pairs)

def list_dot_bracket(db):
    open_stack = []
    out = np.zeros(len(db))
    for i, c in enumerate(db):
        if c == '(':
            open_stack.append(i)
        elif c == ')':
            j = open_stack.pop()
            out[i] = j
            out[j] = i
        elif c == '.':
            out[i] = -1


    return(out)

# Parse an rdat file into a dict with sequence and reactivities
def parse_rdat(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    data_ids = []
    data = {}

    for l in lines:
        if not l.startswith('DATA'):
            continue
        if l.startswith('DATA_ANNOTATION:'):
            if 'datatype:REACTIVITY' in l:
                l = l.split()
                id = l[0].strip('DATA_ANNOTATION:')
                data_ids.append(id)
                data[id] = {'seq' : l[1].strip('sequence:')}
        elif l.startswith('DATA:'):
            l = l.split()
            id = l[0].strip('DATA:')
            if id in data_ids:
                data[id]['react'] = [float(r) for r in l[1:]]

    # Renumber the data keys to match sequence lengths
    # This only works because the original keys are strings and new ones are ints.
    keys = list(data.keys()) # Can't del stuff in dict in a for loop over dictkeys
    for k in keys:
        data[len(data[k]['seq'])] = data[k]
        del(data[k])

    return(data)

# Display in a Forna iframe
def forna_display(seq, struct, cols={}):
    col_str = '\\n'.join([f'{k}:{cols[k]}' for k in cols.keys()])
    # Can we run FORNA in a notebook (I feel like I've done this before...)
    forna_src = f'http://rna.tbi.univie.ac.at/forna/forna.html?id=url/name&sequence={seq}&structure={struct}&colors={col_str}'
    #generate a unique id for our iframe
    return IFrame(forna_src, 1000, 500)
    # That was easy.