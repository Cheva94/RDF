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
from numpy import sqrt, array, zeros, inner, pi, arccos, mean, std

def main():
    start = time()

    name = args.input_file
    center = args.central_atom
    neigh1 = args.neighboring_atom[0]
    neigh2 = args.neighboring_atom[1]
    Rcut = args.Rcut
    Rcut2 = Rcut * Rcut
    convfact = 180 / pi

    print(f'Running bond angle for {center} surrounded by {neigh1} and {neigh2}.')

    xyz = read_csv(name, header = None, delim_whitespace = True,
                    names=['idAt', 'rx', 'ry', 'rz'])
    xyz = xyz.iloc[1:,:].reset_index(drop=True)

    aux = []
    angs = []
    if neigh1 == neigh1 == center:
        nAt = xyz.iloc[:,0].value_counts()[neigh1]
        xyz = xyz[xyz['idAt'] == neigh1].to_numpy()
        r = xyz[:, 1:]

        for i in range(nAt):
            for j in range(i+1, nAt):
                r2 = r[i] - r[j]
                d2 = inner(r2, r2)
                if d2 <= Rcut2:
                    aux.append((i,j,sqrt(d2), r2))

        L = len(aux)
        for m in range(L):
            for n in range(m+1, L):
                if aux[m][0] == aux[n][0]:
                    angulo = convfact * arccos(inner(aux[m][3], aux[n][3]) / (aux[m][2] * aux[n][2]))
                    angs.append(angulo)
                elif aux[m][1] == aux[n][1]:
                    angulo = convfact * arccos(inner(aux[m][3], aux[n][3]) / (aux[m][2] * aux[n][2]))
                    angs.append(angulo)

        A = array(angs)
        print(f'Job done in {(time() - start):.3f} seconds!')

        print(f'The list of bond angles for {neigh1}-{center}-{neigh2} is (in degrees): \n\t{angs} \nThe average bond angle for {neigh1}-{center}-{neigh2} is {mean(A):.2f} degrees with a standard deviation of {std(A):.2f} degrees.')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('input_file', help = "Path to the xyz input file.")

    parser.add_argument('Rcut', type = float, help = "Maximum distance to be \
                        considered.")

    parser.add_argument('central_atom')

    parser.add_argument('neighboring_atom', nargs = 2)

    args = parser.parse_args()

    main()
