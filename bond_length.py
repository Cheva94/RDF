#!/usr/bin/python3.9

'''
    Calculation: bond length.
    Description: Determines the bond length of a pair of atoms for a system when
                given a xsf. Everything is in Angstrom. It takes into account PBC.
    Written by: Ignacio J. Chevallier-Boutell.
    Dated: October, 2021.
'''

import argparse
from pandas import read_csv
from time import time
from numpy import sqrt, inner, mean, std

def PBC(dist, length):
    '''
    Correct a distance using Periodic Boundary Conditions (PBC).
    '''

    return dist - length * int(2*dist/length)

def main():
    start = time()

    name = args.input_file
    at1 = args.atoms[0]
    at2 = args.atoms[1]
    Rcut = args.Rcut
    Rcut2 = Rcut * Rcut
    Rmin = args.Rmin
    if Rmin == None:
        Rmin = 0
    Rmin2 = Rmin * Rmin

    xsf = read_csv(name, header = None, delim_whitespace = True,
                    names=['idAt', 'rx', 'ry', 'rz', 'fx', 'fy', 'fz'])
    Lx, Ly, Lz = float(xsf.iloc[3,0]), xsf.iloc[4,1], xsf.iloc[5,2]
    xyz = xsf.iloc[6:,0:4].reset_index(drop=True)

    ID = []
    POS = []
    BL = []
    if at1 != at2:
        nAt1 = xyz.iloc[:,0].value_counts()[at1]
        nAt2 = xyz.iloc[:,0].value_counts()[at2]

        xyz1 = xyz[xyz['idAt'] == at1]
        xyz2 = xyz[xyz['idAt'] == at2]

        rx1 = xyz1.iloc[:, 1]
        ry1 = xyz1.iloc[:, 2]
        rz1 = xyz1.iloc[:, 3]

        rx2 = xyz2.iloc[:, 1]
        ry2 = xyz2.iloc[:, 2]
        rz2 = xyz2.iloc[:, 3]

        for i in range(nAt1):
            for j in range(nAt2):
                dx = rx1.iloc[i] - rx2.iloc[j]
                dy = ry1.iloc[i] - ry2.iloc[j]
                dz = rz1.iloc[i] - rz2.iloc[j]

                dx = PBC(dx, Lx)
                dy = PBC(dy, Ly)
                dz = PBC(dz, Lz)

                d2 = dx**2 + dy**2 + dz**2
                if ((d2>= Rmin2) and (d2<= Rcut2)):
                    ID.append(f'{rx1.index[i]} - {rx2.index[j]}')
                    POS.append(f'({rx1.iloc[i]}; {ry1.iloc[i]}; {rz1.iloc[i]}) - ({rx2.iloc[j]}; {ry2.iloc[j]}; {rz2.iloc[j]})')
                    BL.append(sqrt(d2))

    else:
        nAt = xyz.iloc[:,0].value_counts()[at1]
        xyz = xyz[xyz['idAt'] == at1]
        id = xyz.iloc[:, 0]
        rx = xyz.iloc[:, 1]
        ry = xyz.iloc[:, 2]
        rz = xyz.iloc[:, 3]

        for i in range(nAt):
            for j in range(i+1, nAt):
                dx = rx.iloc[i] - rx.iloc[j]
                dy = ry.iloc[i] - ry.iloc[j]
                dz = rz.iloc[i] - rz.iloc[j]

                dx = PBC(dx, Lx)
                dy = PBC(dy, Ly)
                dz = PBC(dz, Lz)

                d2 = dx**2 + dy**2 + dz**2
                if ((d2>= Rmin2) and (d2<= Rcut2)):
                    ID.append(f'{rx.index[i]} - {rx.index[j]}')
                    POS.append(f'({rx.iloc[i]}; {ry.iloc[i]}; {rz.iloc[i]}) - ({rx.iloc[j]}; {ry.iloc[j]}; {rz.iloc[j]})')
                    BL.append(sqrt(d2))

    Summary = f'Summary >> {at1}-{at2} = ({mean(BL):.4f} +- {std(BL):.4f}) A ; Count = {len(BL)} <<'

    with open(f'BL_{at1}-{at2}.csv', 'w') as f:
        f.write('==== Bond Length in Angstroms ==== \n\n')
        f.write('Atoms ID, Atoms Position, Distance \n')
        for i in range(len(BL)):
            f.write(f'{ID[i]}, {POS[i]}, {BL[i]:.4f} \n')
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

    parser.add_argument('--Rmin', type = float,
                        help = "Minimum distance to be considered as bond length.")

    args = parser.parse_args()

    main()
