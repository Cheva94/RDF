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

        L = []
        for i in range(nAt):
            for j in range(i+1, nAt):
                r2 = r[i] - r[j]
                d2 = inner(r2, r2)
                if d2 <= Rcut2:
                    L.append((i,j,sqrt(d2), r2))

        angs = []
        for m in range(len(L)):
            for n in range(m+1, len(L)):
                if L[m][0] == L[n][0]:
                    # print(f'match en el primer elemento: {L[m]} con {L[n]}')
                    angulo = convfact * arccos(inner(L[m][3], L[n][3]) / (L[m][2] * L[n][2]))
                    angs.append(angulo)
                elif L[m][1] == L[n][1]:
                    angulo = convfact * arccos(inner(L[m][3], L[n][3]) / (L[m][2] * L[n][2]))
                    angs.append(angulo)

        print(angs)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('input_file', help = "Path to the xyz input file.")

    parser.add_argument('Rcut', type = float, help = "Maximum distance to be \
                        considered.")

    parser.add_argument('c_at', help = "Central atom.")

    parser.add_argument('n_at', nargs = 2, help = "Neighboring atom.")

    args = parser.parse_args()

    main()
