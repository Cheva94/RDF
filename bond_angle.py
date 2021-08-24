#!/usr/local/bin/python3.9

'''
    Calculation: Internal coordniates.
    Description: Determines internal coordinates of a system when given a xyz
                file: bond length, bond angle and . Everything is in Angstrom.
    Written by: Ignacio J. Chevallier-Boutell.
    Dated: August, 2021.
'''

import argparse
from pandas import read_csv
from time import time
from numpy import sqrt, array, zeros, inner#, pi
from os.path import isfile

def main():
    start = time()

    name = args.input_file
    print(f'Running internal coordinates for {name} file.')

    xyz = read_csv(name, header = None, delim_whitespace = True,
                    names=['idAt', 'rx', 'ry', 'rz'])

    idAt = xyz['idAt'].drop_duplicates().to_numpy()
    xyz = xyz.to_numpy()
    nAt = xyz.shape[0]

    dict = {}
    At = xyz[:,0]
    r = xyz[:, 1:]
    for i in range(nAt):
        for j in range(i+1, nAt):
            for k in range(j+1, nAt):
                r12 = r[i] - r[j]
                d12 = inner(r12, r12)
                r23 = r[j] - r[k]
                d23 = inner(r23, r23)

                if d12 <= 4 and d23 <=4:
                    coseno = inner(r12, r23) / sqrt(d12 * d23)
                    dict.setdefault(f'{At[i]}-{At[j]}-{At[k]}', []).append(coseno)

    print(dict.keys())
    # name = name.split('.xyz')[0].split('/')[-1]
    # for at1 in idAt:
    #     for at2 in idAt:
    #         if at1 == at2:
    #             H = zeros(301, dtype=int)
    #             for data in dict.get(f'{at1}-{at2}'):
    #                 H[int(data/0.01)] += 1
    #
    #             with open(f'{name}-bond_length-{at1}-{at2}.csv', 'w') as f:
    #                 for binIdx in range(301):
    #                     if H[binIdx] != 0:
    #                         d = (binIdx + 0.5) * 0.01
    #                         f.write(f'{d:.3f}, {H[binIdx]} \n')
    #         elif at1 != at2:
    #             if isfile(f'{name}-bond_length-{at2}-{at1}.csv'):
    #                 continue
    #             else:
    #                 H = zeros(301, dtype=int)
    #                 for data in dict.get(f'{at1}-{at2}'):
    #                     H[int(data/0.01)] += 1
    #
    #                 for data in dict.get(f'{at2}-{at1}'):
    #                     H[int(data/0.01)] += 1
    #
    #                 with open(f'{name}-bond_length-{at1}-{at2}.csv', 'w') as f:
    #                     for binIdx in range(301):
    #                         if H[binIdx] != 0:
    #                             d = (binIdx + 0.5) * 0.01
    #                             f.write(f'{d:.3f}, {H[binIdx]} \n')

    print(f'Job done in {(time() - start):.3f} seconds!')
    print(f'Output files: {name}-bond_length-ATOMS.csv')
    # print(f'There are {nAtSlab/frames_count:.2f} {at} atoms on average within this slab.')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('input_file', help = "Path to the xyz input file.")
    #
    # parser.add_argument('dh', type = float, help = "Increment to be considered.")
    #
    # parser.add_argument('at', help = "Atom to be analyzed.")
    #
    parser.add_argument('-bL', '--bond_length', action = 'store_true',
                        help = "Calculates bond lengths.")
    #
    # parser.add_argument('-o', '--output_file', help = "Path to the output file. \
    #                     If not given, the default name will be used.")
    #
    # parser.add_argument('-H', '--Hcut', type = float, nargs = 2, default = [0, -1],
    #                     help = "Minimum and maximum heights to be considered.")
    #
    args = parser.parse_args()

    main()
