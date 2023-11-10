import RNA
import numpy as np
from utils import dict_dot_bracket

# According to https://pubs.acs.org/doi/full/10.1021/jacs.0c03105
# dG for each BM step is between 7.4 and 8.9 KbT (4.5-5.5 kcal/mol @ 37)

#penalize opening a base pair by a constant amount
def penalize_barriers(i, j, k, l, d, arg_dict):
    penalty = arg_dict['penalty'] # 20 is the best constant value for this.
    if d in [RNA.DECOMP_PAIR_IL, RNA.DECOMP_PAIR_HP, RNA.DECOMP_PAIR_ML]:
        ref = arg_dict['last']
        if i in ref.keys():
            if ref[i] == j:
                return 0
            else:
                return penalty
            
    return 0

# eval_move is really expensive, so pre-assemble penalties for penalize_barriers_seq
def get_penalties(last, fc):
    last_dict = dict_dot_bracket(last)
    out = {}
    for i, c in enumerate(last):
        if c == '(' or c == ')':
            p = min((i+1, last_dict[i+1]))
            q = max((i+1, last_dict[i+1]))
            out[i+1] = int(fc.eval_move(last, -p, -q) * 100) # eval_move returns kcal/mol, penalties must be in dcal/mol

    return out


#penalize opening a base pair by a percentage of the FE loss
def penalize_barriers_seq(i, j, k, l, d, arg_dict):
    if d in [RNA.DECOMP_PAIR_IL, RNA.DECOMP_PAIR_HP, RNA.DECOMP_PAIR_ML]:
        penalty_percent = arg_dict['penalty_percent']
        ref_dict = arg_dict['last_dict']
        penalty_dict = arg_dict['penalty_dict']

        if i in ref_dict.keys():
            if ref_dict[i] == j:
                return 0
            else:
                diff = penalty_dict[i]
                penalty = int(diff * penalty_percent)
                return penalty
            
    return 0 

# Get the frequency with which a particular base pair is formed
# Optionally also compute the cost of breaking that pair
def pairing_frequency(ensemble, last_fc=None):
    freqs = np.zeros((len(ensemble[0].structure), len(ensemble[0].structure)))
    energies = np.array([e.energy for e in ensemble])
    energies *= -1 * 1.624 # convert to negative KbTs
    probs = np.exp(energies) / np.sum(np.exp(energies))
    for e, prob in zip(ensemble, probs):
        structure = e.structure
        db = dict_dot_bracket(structure)
        for k in db.keys():
            if k > db[k]: # calculating eval_move on i > j returns bad stuff
                continue
            if last_fc == None:
                # Probability that i and j are paired
                freqs[k-1][db[k]-1] += prob # Vienna is 1-indexed, so dict_db is as well
            else:
                # Probability that i and j are paired, multiplied by the energy penalty for breaking them
                move_cost = last_fc.eval_move(structure, -1*k, -1*db[k])
                freqs[k-1][db[k]-1] += prob * move_cost * 100 # eval_move returns kcal/mol, need dcal/mol
                freqs[db[k]-1][k-1] += prob * move_cost * 100 # do it symmetrically so there's a penalty on opening closing nucleotides

    # What we have is:  freqs[i][j] = p(i_paired_to_j)
    # What we want is:  freqs[i][j] = p(i_paired_to_not_j)
    freqs -= np.sum(freqs, axis=1)[:, np.newaxis] # negative penalties can happen!
    freqs *= -1

    return freqs

# Penalty based on the frequency with which i,j are not paired in the ensemble
def penalize_barriers_ensemble(i, j, k, l, d, arg_dict):
    if d in [RNA.DECOMP_PAIR_IL, RNA.DECOMP_PAIR_HP, RNA.DECOMP_PAIR_ML]:
        freqs = arg_dict['freqs']
        penalty = arg_dict['penalty']
        return int(freqs[i-1][j-1] * penalty)

    return 0

# All the constraint funcitons have the same arguments so I can call them en-mass
def no_constraint(seq, _, _2, md=RNA.md()):
    fc = RNA.fold_compound(seq, md)
    return fc.mfe()

def shape_constraint(seq, reactivities, md=RNA.md()):
    fc = RNA.fold_compound(seq, md)
    #fc.sc_add_SHAPE_deigan(reactivities, 2.6, -0.8)
    fc.sc_add_SHAPE_zarringhalam(reactivities, 0.8, 0.5, 'M')
    return fc.mfe()

def constant_penalty(seq, penalty, last_structure, md=RNA.md()):
    fc = RNA.fold_compound(seq, md)
    step_info = {
            'last' : dict_dot_bracket(last_structure),
            'penalty' : int(penalty)
        }
    fc.sc_add_f(penalize_barriers)
    fc.sc_add_data(step_info)

    return fc.mfe()

def sequence_dependent_penalty(seq, penalty_percent, last_structure, md=RNA.md()):
    fc = RNA.fold_compound(seq, md)
    if last_structure != '':
        fc_last = RNA.fold_compound(seq[:len(last_structure)], md)
        penalty_dict = get_penalties(last_structure, fc_last)
        #p_list.extend([v for v in penalty_dict.values()]) #was used to get average penalty
    else:
        penalty_dict = {}
    step_info = {
        'last_dict' : dict_dot_bracket(last_structure),
        'penalty_dict' : penalty_dict,
        'penalty_percent' : penalty_percent
        }
    fc.sc_add_f(penalize_barriers_seq)
    fc.sc_add_data(step_info)

    return fc.mfe()

# This didn't work very well and made the function swapping complicated with its extra argument
#def hierarchical_fold(seq, penalty, last_structure, span):
#    mds = RNA.md()
#    mds.temperature = 37
#    mds.max_bp_span = span
#    step_info = {
#        'last' : dict_dot_bracket(last_structure),
#        'penalty' : int(penalty)
#    }
#    fc = RNA.fold_compound(seq, mds)
#    fc.sc_add_f(penalize_barriers)
#    fc.sc_add_data(step_info)
#
#    return fc.mfe()

def constant_ensemble_penalty(seq, penalty, last_ensemble, md=RNA.md()):
    if last_ensemble != []:
        last_freqs = pairing_frequency(last_ensemble)
    else:
        last_freqs = np.zeros(len(seq))

    freqs = np.zeros((len(seq), len(seq)))
    freqs[:len(last_freqs), :len(last_freqs)] = last_freqs
        
    step_info = {
        'freqs' : freqs,
        'penalty' : int(penalty)
    }
    fc = RNA.fold_compound(seq, md)
    fc.sc_add_f(penalize_barriers_ensemble)
    fc.sc_add_data(step_info)
    return fc.subopt(500)

def sequence_dependent_ensemble_penalty(seq, percent, last_ensemble, md=RNA.md()):
    if last_ensemble != []:
        fc_last = RNA.fold_compound(seq[:len(last_ensemble[0].structure)], md)
        last_freqs = pairing_frequency(last_ensemble, fc_last)
    else:
        last_freqs = np.zeros(len(seq))

    freqs = np.zeros((len(seq), len(seq)))
    freqs[:len(last_freqs), :len(last_freqs)] = last_freqs
        
    step_info = {
        'freqs' : freqs,
        'penalty' : percent
    }
    fc = RNA.fold_compound(seq, md)
    fc.sc_add_f(penalize_barriers_ensemble)
    fc.sc_add_data(step_info)
    return fc.subopt(500)