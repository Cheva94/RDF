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
from numpy import sqrt, inner, mean, std

def main():
    start = time()

    name = args.input_file
    at1 = args.atoms[0]
    at2 = args.atoms[1]
    Rcut = args.Rcut
    Rcut2 = Rcut * Rcut

    xyz = read_csv(name, header = None, delim_whitespace = True,
                    names=['idAt', 'rx', 'ry', 'rz'])
    xyz = xyz.iloc[1:,:].reset_index(drop=True)

    ID = []
    BL = []
    if at1 != at2:
        nAt1 = xyz.iloc[:,0].value_counts()[at1]
        nAt2 = xyz.iloc[:,0].value_counts()[at2]

        xyz1 = xyz[xyz['idAt'] == at1].to_numpy()
        xyz2 = xyz[xyz['idAt'] == at2].to_numpy()

        id1 = xyz1[:, 0]
        r1 = xyz1[:, 1:]
        id2 = xyz2[:, 0]
        r2 = xyz2[:, 1:]

        for i in range(nAt1):
            for j in range(nAt2):
                d2 = r1[i] - r2[j]
                d2 = inner(d2, d2)
                if d2 <= Rcut2:
                    ID.append(f'{i} - {j}')
                    BL.append(sqrt(d2))

    else:
        nAt = xyz.iloc[:,0].value_counts()[at1]
        xyz = xyz[xyz['idAt'] == at1].to_numpy()
        id = xyz[:, 0]
        r = xyz[:, 1:]

        for i in range(nAt):
            for j in range(i+1, nAt):
                d2 = r[i] - r[j]
                d2 = inner(d2, d2)
                if d2 <= Rcut2:
                    ID.append(f'{i} - {j}')
                    BL.append(sqrt(d2))

    Summary = f'Summary >> {at1}-{at2} = ({mean(BL):.4f} +- {std(BL):.4f}) A ; Count = {len(BL)} <<'

    with open(f'{at1}-{at2}.csv', 'w') as f:
        f.write('==== Bond Length in Angstroms ==== \n\n')
        f.write('Atoms ID, Distance \n')
        for i in range(len(BL)):
            f.write(f'{ID[i]}, {BL[i]:.4f} \n')
        # f.write(f'\n Summary: {at1}-{at2} = ({mean(BL):.4f} +- {std(BL):.4f}) A ; Count = {len(BL)}')
        f.write(f'\n {Summary}')

    print(f'Job done in {(time() - start):.3f} seconds!')
    print(f'Output file: {at1}-{at2}.csv')
    print(f'{Summary}')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('input_file', help = "Path to the xyz input file.")

    parser.add_argument('Rcut', type = float, help = "Maximum distance to be \
                        considered as bond length.")

    parser.add_argument('atoms', nargs = 2, help = "Atoms to be analyzed.")

    args = parser.parse_args()

    main()
