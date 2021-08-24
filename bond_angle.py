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
    aux2 = []
    angs = []

    if neigh1 == neigh2 == center:
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
                if aux[m][0] == aux[n][0] or aux[m][1] == aux[n][1]:
                    angulo = convfact * arccos(inner(aux[m][3], aux[n][3]) / (aux[m][2] * aux[n][2]))
                    angs.append(angulo)

        A = array(angs)
        print(f'Job done in {(time() - start):.3f} seconds!')
        print(f'The list of bond angles for {neigh1}-{center}-{neigh2} is (in degrees): \n\t{angs} \nThe average bond angle for {neigh1}-{center}-{neigh2} is {mean(A):.1f} degrees with a standard deviation of {std(A):.1f} degrees.')

    if neigh1 == neigh2 != center:
        nCenter = xyz.iloc[:,0].value_counts()[center]
        nNeigh = xyz.iloc[:,0].value_counts()[neigh1]

        xyzCenter = xyz[xyz['idAt'] == center].to_numpy()
        xyzNeigh = xyz[xyz['idAt'] == neigh1].to_numpy()

        rCenter = xyzCenter[:, 1:]
        rNeigh = xyzNeigh[:, 1:]

        for i in range(nCenter):
            for j in range(nNeigh):
                r2 = rNeigh[j] - rCenter[i]
                d2 = inner(r2, r2)
                if d2 <= Rcut2:
                    aux.append((i,j,sqrt(d2), r2))

        L = len(aux)
        for m in range(L):
            for n in range(m+1, L):
                if aux[m][0] == aux[n][0]:
                    angulo = convfact * arccos(inner(aux[m][3], aux[n][3]) / (aux[m][2] * aux[n][2]))
                    angs.append(angulo)

        A = array(angs)
        print(f'Job done in {(time() - start):.3f} seconds!')
        print(f'The list of bond angles for {neigh1}-{center}-{neigh2} is (in degrees): \n\t{angs} \nThe average bond angle for {neigh1}-{center}-{neigh2} is {mean(A):.1f} degrees with a standard deviation of {std(A):.1f} degrees.')

    if neigh1 == center != neigh2:
        nCenter = xyz.iloc[:,0].value_counts()[center]
        nNeigh = xyz.iloc[:,0].value_counts()[neigh2]

        xyzCenter = xyz[xyz['idAt'] == center].to_numpy()
        xyzNeigh = xyz[xyz['idAt'] == neigh2].to_numpy()

        rCenter = xyzCenter[:, 1:]
        rNeigh = xyzNeigh[:, 1:]

        for i in range(nCenter):
            for j in range(nNeigh):
                r2 = rNeigh[j] - rCenter[i]
                d2 = inner(r2, r2)
                if d2 <= Rcut2:
                    aux.append((i,j,sqrt(d2), r2))

        for i in range(nCenter):
            for j in range(i+1, nCenter):
                r2 = rCenter[j] - rCenter[i]
                d2 = inner(r2, r2)
                if d2 <= Rcut2:
                    aux2.append((i,j,sqrt(d2), r2))

        L = len(aux)
        L2 = len(aux2)
        for m in range(L):
            for n in range(L2):
                if aux[m][0] == aux2[n][0]:
                    angulo = convfact * arccos(inner(aux[m][3], aux2[n][3]) / (aux[m][2] * aux2[n][2]))
                    angs.append(angulo)
                elif aux[m][0] == aux2[n][1]:
                    angulo = convfact * arccos(inner(aux[m][3], -aux2[n][3]) / (aux[m][2] * aux2[n][2]))
                    angs.append(angulo)

        A = array(angs)
        print(f'Job done in {(time() - start):.3f} seconds!')
        print(f'The list of bond angles for {neigh1}-{center}-{neigh2} is (in degrees): \n\t{angs} \nThe average bond angle for {neigh1}-{center}-{neigh2} is {mean(A):.1f} degrees with a standard deviation of {std(A):.1f} degrees.')

    if neigh2 == center != neigh1:
        nCenter = xyz.iloc[:,0].value_counts()[center]
        nNeigh = xyz.iloc[:,0].value_counts()[neigh1]

        xyzCenter = xyz[xyz['idAt'] == center].to_numpy()
        xyzNeigh = xyz[xyz['idAt'] == neigh1].to_numpy()

        rCenter = xyzCenter[:, 1:]
        rNeigh = xyzNeigh[:, 1:]

        for i in range(nCenter):
            for j in range(nNeigh):
                r2 = rNeigh[j] - rCenter[i]
                d2 = inner(r2, r2)
                if d2 <= Rcut2:
                    aux.append((i,j,sqrt(d2), r2))

        for i in range(nCenter):
            for j in range(i+1, nCenter):
                r2 = rCenter[j] - rCenter[i]
                d2 = inner(r2, r2)
                if d2 <= Rcut2:
                    aux2.append((i,j,sqrt(d2), r2))

        L = len(aux)
        L2 = len(aux2)
        for m in range(L):
            for n in range(L2):
                if aux[m][0] == aux2[n][0]:
                    angulo = convfact * arccos(inner(aux[m][3], aux2[n][3]) / (aux[m][2] * aux2[n][2]))
                    angs.append(angulo)
                elif aux[m][0] == aux2[n][1]:
                    angulo = convfact * arccos(inner(aux[m][3], -aux2[n][3]) / (aux[m][2] * aux2[n][2]))
                    angs.append(angulo)

        A = array(angs)
        print(f'Job done in {(time() - start):.3f} seconds!')
        print(f'The list of bond angles for {neigh1}-{center}-{neigh2} is (in degrees): \n\t{angs} \nThe average bond angle for {neigh1}-{center}-{neigh2} is {mean(A):.1f} degrees with a standard deviation of {std(A):.1f} degrees.')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('input_file', help = "Path to the xyz input file.")

    parser.add_argument('Rcut', type = float, help = "Maximum distance to be \
                        considered.")

    parser.add_argument('central_atom')

    parser.add_argument('neighboring_atom', nargs = 2)

    args = parser.parse_args()

    main()
