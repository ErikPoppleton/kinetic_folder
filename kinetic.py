import RNA
import os
import numpy as np
from scipy.spatial import distance
import matplotlib.pyplot as plt
import argparse
from utils import dict_dot_bracket

# If something in the full seq is paired to something which doesn't exist in the subseq
# Then I don't think it counts as a missfold if that nucleotide is unpaired in the subseq fold
def partial_ideal(ref):
    # First we identify unclosed opening brackets
    ref_transformed = []
    pstack = []
    for i, c in enumerate(ref):
        ref_transformed.append(c)
        if c == '(':
            pstack.append(i)
        if c == ')':
            pstack.pop()

    # And pretend those are unpaired.
    for i in pstack:
        ref_transformed[i] = '.'
    
    return(ref_transformed)



# Energy correction which penalizes base pairs which were not present in the last iteration
# The penalty is equal to ...% of the dinucleotide step which needs to be broken
# Cody suggests 1.5-2.5 kcal
def penalize_barriers(i, j, k, l, d, arg_dict):
    penalty = 500 # WHAT SHOULD THIS BE???
    if d in [RNA.DECOMP_PAIR_IL, RNA.DECOMP_PAIR_HP, RNA.DECOMP_PAIR_ML]:
        ref = arg_dict['last']
        if i in ref.keys():
            if ref[i] == j:
                return 0
            else:
                return penalty
            
    return 0
        

files = os.listdir('.')
files = [f for f in files if 'design' in f]

md = RNA.md()
md.temperature = 37

for f in files:
    with open(f, 'r') as design:
        name = design.readline().strip()
        target = design.readline().strip()
        seq = design.readline().strip()

    # No pseudoknots in RNAsubopt
    target = target.replace('[', '.')
    target = target.replace(']', '.')

    # Create subsequences, each 10 nt longer than the last.  Fold.
    # Then compare the distance between the fold of the subseq to the truncated fold of the whole structure.
    # This tries both a fresh fold of the whole subsequence, as well as a fold with soft constraints based on the MFE of the last iteration.
    step_size = 10
    boundspace = np.arange(10, len(seq)+step_size, step_size)

    mins_comp = np.empty_like(boundspace)
    maxes_comp = np.empty_like(boundspace)
    meds_comp = np.empty_like(boundspace)
    mfes_comp = np.empty_like(boundspace)
    n_sub_comp = np.empty_like(boundspace)
    mins_con = np.empty_like(boundspace)
    maxes_con = np.empty_like(boundspace)
    meds_con = np.empty_like(boundspace)
    mfes_con = np.empty_like(boundspace)
    n_sub_con = np.empty_like(boundspace)
    mfe_diff = np.empty_like(boundspace)
    last_structure = {}

    for i, bound in enumerate(boundspace):
        subseq = seq[:bound]
        subtarget = target[:bound]
        ref = partial_ideal(subtarget)

        ###################################
        #          complete fold          #
        ###################################
        fc = RNA.fold_compound(subseq, md)
        competing = fc.subopt(100)
        
        energies = [c.energy for c in competing]
        structures = [c.structure for c in competing]
        distances = [int(np.round(distance.hamming(list(s), ref) * len(ref))) for s in structures]
        min_idx = np.argmin(energies)
        mfe_comp = structures[min_idx]

        mins_comp[i] = min(distances)
        maxes_comp[i] = max(distances)
        meds_comp[i] = np.median(distances)
        mfes_comp[i] = distances[min_idx]
        n_sub_comp[i] = len(competing)

        ###################################
        #         constrained fold        #
        ###################################
        fc = RNA.fold_compound(subseq, md)

        step_info = {
            'last' : last_structure
        }

        fc.sc_add_f(penalize_barriers)
        fc.sc_add_data(step_info)

        competing = fc.subopt(100)

        distances = [int(np.round(distance.hamming(list(c.structure), ref) * len(ref))) for c in competing]
        energies = [c.energy for c in competing]
        structures = [c.structure for c in competing]

        min_idx = np.argmin(energies)
        mfe_con = structures[min_idx]
        last_structure = dict_dot_bracket(mfe_con)

        mins_con[i] = min(distances)
        maxes_con[i] = max(distances)
        meds_con[i] = np.median(distances)
        mfes_con[i] = distances[min_idx]
        n_sub_con[i] = len(competing)


        mfe_diff[i] = int(np.round(distance.hamming(list(mfe_comp), list(mfe_con)) * len(mfe_con)))
        print(f"Curr_length: {len(subseq)}, Comp_subs: {n_sub_comp[i]}, Con_subs: {n_sub_con[i]}, MFE_dist: {mfe_diff[i]}")


    # let's plot the min, med and max and see that path to the final structure.
    fig, ax = plt.subplots(1, 2, figsize=(10, 5))
    a02 = ax[0].twinx()
    ln1 = a02.plot(boundspace, n_sub_comp, 'r--', label='n structures')
    ln2 = ax[0].plot(boundspace, maxes_comp, label='max dist', c='deepskyblue')
    ln3 = ax[0].plot(boundspace, meds_comp, label='median dist', c='royalblue')
    ln4 = ax[0].plot(boundspace, mfes_comp, label='mfe dist', c='g')
    ln5 = ax[0].plot(boundspace, mins_comp, label='min dist', c='darkblue')
    ax[0].set_title("Complete fold")
    ax[0].set_xlabel("subsequence length")

    lines = ln1 + ln2 + ln3 + ln4 + ln5
    labs = [l.get_label() for l in lines]
    ax[0].legend(lines, labs)

    a12 = ax[1].twinx()
    ln1 = a12.plot(boundspace, n_sub_con, 'r--', label='n structures')
    ln2 = ax[1].plot(boundspace, maxes_con, label='max dist', c='deepskyblue')
    ln3 = ax[1].plot(boundspace, meds_con, label='median dist', c='royalblue')
    ln4 = ax[1].plot(boundspace, mfes_con, label='mfe dist', c='g')
    ln5 = ax[1].plot(boundspace, mins_con, label='min dist', c='darkblue')
    ax[1].set_xlabel("Subsequence length")
    ax[1].set_title("Constrained fold")

    #lines = ln1 + ln2 + ln3 + ln4 + ln5
    #labs = [l.get_label() for l in lines]
    #ax[1].legend(lines, labs)

    # make first two plots share yscale
    y_min = min([ax[0].get_ylim()[0], ax[1].get_ylim()[0]])
    y_max = max([ax[0].get_ylim()[1], ax[1].get_ylim()[1]])
    ax[0].set_ylim((y_min, y_max))
    ax[1].set_ylim((y_min, y_max))
    ax[0].set_ylabel("Hamming distance")
    ax[0].spines[['right', 'top']].set_visible(False)
    a02.spines[['right', 'top']].set_visible(False)
    a02.set_yticks([])
    a12.set_ylabel("Number of substructures")
    ax[1].spines[['left', 'top']].set_visible(False)
    a12.spines[['left', 'top']].set_visible(False)
    ax[1].set_yticks([])

    plt.tight_layout()
    plt.savefig(f.strip('design.txt')+'kinetic.png', dpi = 200)

    fig, ax = plt.subplots()
    ax.plot(boundspace, mfe_diff)
    ax.set_ylabel("Hamming distance")
    ax.set_xlabel("Subsequence length")
    ax.set_title("Difference in MFE structure")
    plt.tight_layout()
    plt.savefig(f.strip('design.txt')+'fold_diff.png', dpi = 200)
