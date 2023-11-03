import RNA
import os
import numpy as np
from scipy.spatial import distance
import matplotlib.pyplot as plt
import argparse

KCAL_TO_KBT = 1.624 # 1kcal/mol = 1.624 kBT @310.15K.

files = os.listdir('.')
files = [f for f in files if 'design' in f]

md = RNA.md()
md.temperature = 37

parser = argparse.ArgumentParser()
parser.add_argument('delta', type=float, help='Window for suboptimal structure prediction (in kBT)')
args = parser.parse_args()

delta = args.delta #1kBT = 1.023e-21 cal @310.15K.  Vienna is in kcal/mol
delta = delta * KCAL_TO_KBT # 1kcal/mol = 1.624 kBT @310.15K.
print(f"{args.delta} KbT = {delta:.3} kcal/mol")
delta = delta * 100 # subopt input is in dcal/mol
delta = int(delta) # ViennaRNA input must be an int

#1 kBT -> 2x more likely
#2 kBT -> 7x more likely
#3 kBT -> 20x more likely
#4 kBT -> 55x more likely
#5 kBT -> 150x more likely
#6 kBT -> 403x more likely

for f in files:
    print(f"Working on {f}")
    with open(f, 'r') as design:
        name = design.readline().strip()
        target = design.readline().strip()
        seq = design.readline().strip()

    # No pseudoknots in RNAsubopt
    target = target.replace('[', '.')
    target = target.replace(']', '.')

    fc = RNA.fold_compound(seq, md)

    competing = fc.subopt(delta) 

    print(f"Identified {len(competing)} structures.")

    distances = np.array([int(np.round(distance.hamming(list(target), list(c.structure)) * len(target))) for c in competing])
    energies = np.array([c.energy for c in competing])

    # minimum energy for each distance ended up really noisy because the good structures are evenly spaced 2 apart.
    # get the minimum of each pair and use that.
    distance_values = list(set(distances))
    distance_values.sort()
    min_energies = np.empty_like(distance_values)
    min_energies.dtype = float 
    for i, d in enumerate(distance_values):
        distance_mask = distances==d
        min_energies[i] = min(energies[distance_mask])

    lineX = []
    lineY = []
    for d, m in zip(zip(*[iter(distance_values)]*2), zip(*[iter(min_energies)]*2)): #that's pythonic, apparently.
        min_idx = min(range(len(m)), key=m.__getitem__)
        lineX.append(d[min_idx])
        lineY.append(m[min_idx])


    # Plot all the distance/energy pairs and put a bounding line on the min energies.
    fig, ax = plt.subplots()
    ax.scatter(distances, energies, alpha = 0.1)
    ax.plot(lineX, lineY, 'b--')
    ax.set_ylabel('Free energy (kcal/mol)')
    ax.set_xlabel('Hamming distance')
    ax.set_xticks(np.arange(min(distances), max(distances)+1, 1))
    plt.tight_layout()
    plt.savefig(f.strip('design.txt')+'ensemble.png', dpi = 200)

    ## pick the structure with the furthest distance and perform a tree search looking for the target structure.
    #max_distance = max(distances)
    #distance_mask = distances==max_distance


     