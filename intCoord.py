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
from numpy import sqrt, array, zeros #, pi

def main():
    start = time()

    name = args.input_file

    xyz = read_csv(name, header = None, delim_whitespace = True,
                    names=['idAt', 'rx', 'ry', 'rz'])

    idAt = xyz['idAt'].drop_duplicates().to_numpy()
    xyz = xyz.to_numpy()
    nAt = xyz.shape[0]
    At = array(xyz[:,0])

    rx1 = array(xyz[:,1])
    ry1 = array(xyz[:,2])
    rz1 = array(xyz[:,3])

    rx2 = array(xyz[:,1])
    ry2 = array(xyz[:,2])
    rz2 = array(xyz[:,3])

    dict = {}
    for i in range(nAt):
        for j in range(i+1, nAt):
            dx = rx1[i] - rx2[j]
            dy = ry1[i] - ry2[j]
            dz = rz1[i] - rz2[j]

            d2 = dx**2 + dy**2 + dz**2

            if d2 <= 9:
                dict.setdefault(f'{At[i]}-{At[j]}', []).append(sqrt(d2))

    # output_file = args.output_file
    # if output_file == None:
    # output_file = f"{name.split('.xyz')[0].split('/')[-1]}-bond_length"

    name = name.split('.xyz')[0].split('/')[-1]
    for at1 in idAt:
        for at2 in idAt:
            if at1 == at2:
                H = zeros(301, dtype=int)
                for data in dict.get(f'{at1}-{at2}'):
                    H[int(data/0.01)] += 1

                with open(f'{name}-bond_length-{at1}-{at2}.csv', 'w') as f:
                    for binIdx in range(301):
                        if H[binIdx] != 0:
                            d = (binIdx + 0.5) * 0.01
                            f.write(f'{d:.3f}, {H[binIdx]} \n')
            elif at1 != at2:
                H = zeros(301, dtype=int)
                for data in dict.get(f'{at1}-{at2}'):
                    H[int(data/0.01)] += 1

                with open(f'{name}-bond_length-{at1}-{at2}.csv', 'w') as f:
                    for binIdx in range(301):
                        if H[binIdx] != 0:
                            d = (binIdx + 0.5) * 0.01
                            f.write(f'{d:.3f}, {H[binIdx]} \n')
            #     dict_end.setdefault(f'{at1}-{at2}', []).append(dict.get(f'{at1}-{at2}'))
            #     dict_end.setdefault(f'{at1}-{at2}', []).append(dict.get(f'{at2}-{at1}'))

    # print(f'Job done in {(time() - start):.3f} seconds!')
    # print(f'Output file: {output_file}.csv')
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
