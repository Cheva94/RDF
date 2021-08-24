#!/usr/local/bin/python3.9

'''
    Calculation: bond length.
    Description: Determines the bond length of a pair of atoms for a system when
                given a xyz. Everything is in Angstrom.
    Written by: Ignacio J. Chevallier-Boutell.
    Dated: August, 2021.
'''

import argparse
from pandas import read_csv
from time import time
from numpy import sqrt, array, zeros, inner

def main():
    start = time()

    name = args.input_file
    at1 = args.atoms[0]
    at2 = args.atoms[1]
    Rcut = args.Rcut
    increment = 0.01

    nBin = int(Rcut/increment) + 1
    H = zeros(nBin, dtype = int)
    Rcut2 = Rcut * Rcut

    print(f'Running bond length between {at1} and {at2}.')

    xyz = read_csv(name, header = None, delim_whitespace = True,
                    names=['idAt', 'rx', 'ry', 'rz'])
    xyz = xyz.iloc[1:,:].reset_index(drop=True)

    if at1 != at2:
        nAt1 = xyz.iloc[:,0].value_counts()[at1]
        nAt2 = xyz.iloc[:,0].value_counts()[at2]

        xyz1 = xyz[xyz['idAt'] == at1].to_numpy()
        xyz2 = xyz[xyz['idAt'] == at2].to_numpy()

        r1 = xyz1[:, 1:]
        r2 = xyz2[:, 1:]

        for i in range(nAt1):
            for j in range(nAt2):
                d2 = r1[i] - r2[j]
                d2 = inner(d2, d2)
                if d2 <= Rcut2:
                    binIdx = int(sqrt(d2)/increment)
                    H[binIdx] += 1

        name = name.split('.xyz')[0].split('/')[-1]
        name = f'{name}-bond_length-{at1}-{at2}.csv'
        with open(name, 'w') as f:
            for binIdx in range(nBin):
                if H[binIdx] != 0:
                    d = (binIdx + 0.5) * increment
                    f.write(f'{d:.3f}, {H[binIdx]} \n')

        print(f'Job done in {(time() - start):.3f} seconds!')
        print(f'Output files: {name}')

    else:
        nAt = xyz.iloc[:,0].value_counts()[at1]
        xyz = xyz[xyz['idAt'] == at1].to_numpy()
        r = xyz[:, 1:]

        for i in range(nAt):
            for j in range(i+1, nAt):
                d2 = r[i] - r[j]
                d2 = inner(d2, d2)
                if d2 <= Rcut2:
                    binIdx = int(sqrt(d2)/increment)
                    H[binIdx] += 1

        name = name.split('.xyz')[0].split('/')[-1]
        name = f'{name}-bond_length-{at1}-{at1}.csv'
        with open(name, 'w') as f:
            for binIdx in range(nBin):
                if H[binIdx] != 0:
                    d = (binIdx + 0.5) * increment
                    f.write(f'{d:.3f}, {H[binIdx]} \n')

        print(f'Job done in {(time() - start):.3f} seconds!')
        print(f'Output files: {name}')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('input_file', help = "Path to the xyz input file.")

    parser.add_argument('Rcut', type = float, help = "Maximum distance to be \
                        considered.")

    parser.add_argument('atoms', nargs = 2, help = "Atoms to be analyzed.")

    args = parser.parse_args()

    main()
