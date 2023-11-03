import RNA
from utils import dict_dot_bracket

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
            out[i+1] = int(fc.eval_move(last, -p, -q) * 100)

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