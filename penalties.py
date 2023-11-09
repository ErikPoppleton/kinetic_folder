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

# eval_move is really expensive, so pre-assemble penalties 
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

def pairing_frequency(ensemble, last_fc=None):
    freqs = np.zeros((len(ensemble[0].structure), len(ensemble[0].structure)))
    energies = np.array([e.energy for e in ensemble])
    energies *= -1 * 1.624 # convert to negative KbTs
    probs = np.exp(energies) / np.sum(np.exp(energies))
    for e, prob in zip(ensemble, probs):
        structure = e.structure
        print(structure)
        db = dict_dot_bracket(structure)
        print(db)
        for k in db.keys():
            if last_fc == None:
                # Probability that i and j are paired
                freqs[k-1][db[k]-1] += prob # Vienna is 1-indexed, so dict_db is as well
            else:
                # Probability that i and j are paired, multiplied by the energy penalty for breaking them
                freqs[k-1][db[k]-1] += prob * last_fc.eval_move(structure, -1*k, -1*db[k]) * 100

    # What we have is:  freqs[i][j] = p(i_paired_to_j)
    # What we want is:  freqs[i][j] = p(i_paired_to_not_j)
    freqs -= np.sum(freqs, axis=1)[:, np.newaxis] # something isn't right here. Ending up with negative penalties
    freqs *= -1

    return freqs

def penalize_barriers_ensemble(i, j, k, l, d, arg_dict):
    if d in [RNA.DECOMP_PAIR_IL, RNA.DECOMP_PAIR_HP, RNA.DECOMP_PAIR_ML]:
        freqs = arg_dict['freqs']
        penalty = arg_dict['penalty']
        return int(freqs[i-1][j-1] * penalty)

    return 0