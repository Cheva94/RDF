#!/usr/local/bin/python3.9

'''
    Calculation: bond angle.
    Description: Determines the bond angle of a triple of atoms for a system when
                given a xyz. Everything is in Angstrom.
    Written by: Ignacio J. Chevallier-Boutell.
    Dated: August, 2021.
'''

import argparse
from pandas import read_csv
from time import time
from numpy import sqrt, array, zeros, inner, pi, arccos

def main():
    start = time()

    name = args.input_file
    at2 = args.c_at
    at1 = args.n_at[0]
    at3 = args.n_at[1]
    Rcut = args.Rcut
    increment = 0.1

    nBin = int((190 - 80)/increment) + 1
    H = zeros(nBin, dtype = int)
    Rcut2 = Rcut * Rcut
    convfact = 180 / pi

    print(f'Running bond angle for {at2} surrounded by {at1} and {at3}.')

    xyz = read_csv(name, header = None, delim_whitespace = True,
                    names=['idAt', 'rx', 'ry', 'rz'])
    xyz = xyz.iloc[1:,:].reset_index(drop=True)

    if at1 == at2 == at3:
        nAt = xyz.iloc[:,0].value_counts()[at1]
        xyz = xyz[xyz['idAt'] == at1].to_numpy()
        r = xyz[:, 1:]
        for i in range(nAt):
            for j in range(i+1, nAt):
                for k in range(j+1, nAt):
                    r12 = r[i] - r[j]
                    d12 = inner(r12, r12)
                    r23 = r[k] - r[j]
                    d23 = inner(r23, r23)
                    print(d12)
                    print(d23)
                    if d12 <= Rcut2 and d23 <= Rcut2:
                        angulo = convfact * arccos(inner(r12, r23) / sqrt(d12 * d23))
                        print(angulo)
                        if angulo < 80:
                            continue
                        else:
                            binIdx = int((angulo - 80) / increment)
                            H[binIdx] += 1

        name = name.split('.xyz')[0].split('/')[-1]
        name = f'{name}-bond_angle-{at1}-{at1}-{at1}.csv'
        with open(name, 'w') as f:
            for binIdx in range(nBin):
                if H[binIdx] != 0:
                    ang = (binIdx + 0.5) * increment + 80
                    f.write(f'{ang:.3f}, {H[binIdx]} \n')

        print(f'Job done in {(time() - start):.3f} seconds!')
        print(f'Output files: {name}')


    # nAt1 = xyz.iloc[:,0].value_counts()[at1]
    # nAt2 = xyz.iloc[:,0].value_counts()[at2]
    # nAt3 = xyz.iloc[:,0].value_counts()[at3]
    #
    # xyz1 = xyz[xyz['idAt'] == at1].to_numpy()
    # xyz2 = xyz[xyz['idAt'] == at2].to_numpy()
    # xyz3 = xyz[xyz['idAt'] == at3].to_numpy()
    #
    # r1 = xyz1[:, 1:]
    # r2 = xyz2[:, 1:]
    # r3 = xyz3[:, 1:]

    # nAt = xyz.shape[0]
    # At = xyz[:,0]
    # r = xyz[:, 1:]
    # dict = {}
    # for i in range(nAt):
    #     for j in range(nAt):
    #         for k in range(nAt):
    #             # if At[j] != at2:
    #             #     continue
    #             # elif At[i] == at1 and At[k] == at3:
    #             #     0
    #             # elif At[i] == at3 and At[k] == at3:
    #
    #             r12 = r[i] - r[j]
    #             d12 = inner(r12, r12)
    #             r23 = r[k] - r[j]
    #             d23 = inner(r23, r23)
    #
    #             if d12 <= Rcut2 and d23 <= Rcut2:
    #                 coseno = inner(r12, r23) / sqrt(d12 * d23)
    #                 angulo = arccos(coseno) * 180 / pi
    #                 # dict.setdefault(f'{At[i]}-{At[j]}-{At[k]}', []).append(coseno)
    #                 dict.setdefault(f'{At[i]}-{At[j]}-{At[k]}', []).append(angulo)

    # print(dict)
#
    # print(dict.keys())
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

    # print(f'Job done in {(time() - start):.3f} seconds!')
    # print(f'Output files: {name}-bond_length-ATOMS.csv')
    # # print(f'There are {nAtSlab/frames_count:.2f} {at} atoms on average within this slab.')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('input_file', help = "Path to the xyz input file.")

    parser.add_argument('Rcut', type = float, help = "Maximum distance to be \
                        considered.")

    parser.add_argument('c_at', help = "Central atom.")

    parser.add_argument('n_at', nargs = 2, help = "Neighboring atom.")

    args = parser.parse_args()

    main()
