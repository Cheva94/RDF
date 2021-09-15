#!/usr/local/bin/python3.9

'''
    Calculation: bond angle.
    Description: Determines the bond angle of a triple of atoms for a system when
                given a xyz. Everything is in BAtrom.
    Written by: Ignacio J. Chevallier-Boutell.
    Dated: August, 2021.
'''

import argparse
from pandas import read_csv
from time import time
from numpy import sqrt, inner, pi, arccos, mean, std

def PBC(dist, length):
    '''
    Correct a distance using Periodic Boundary Conditions (PBC).
    '''

    return dist - length * int(2*dist/length)

def main():
    start = time()

    name = args.input_file
    center = args.central_atom
    neigh1 = args.neighboring_atom[0]
    neigh2 = args.neighboring_atom[1]
    Rcut = args.Rcut
    Rcut2 = Rcut * Rcut
    convfact = 180 / pi
    Rmin = args.Rmin
    if Rmin == None:
        Rmin = 0
    Rmin2 = Rmin * Rmin

    xsf = read_csv(name, header = None, delim_whitespace = True,
                    names=['idAt', 'rx', 'ry', 'rz', 'fx', 'fy', 'fz'])
    Lx, Ly, Lz = float(xsf.iloc[3,0]), xsf.iloc[4,1], xsf.iloc[5,2]
    xyz = xsf.iloc[6:,0:4].reset_index(drop=True)

    BA = []
    ID = []
    aux1 = []
    aux2 = []

    if neigh1 == neigh2 == center:
        nAt = xyz.iloc[:,0].value_counts()[center]
        xyz = xyz[xyz['idAt'] == center]
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
                r2 = [dx, dy, dz]

                d2 = dx**2 + dy**2 + dz**2
                if ((d2>= Rmin2) and (d2<= Rcut2)):
                    aux1.append((rx.index[i],rx.index[j],sqrt(d2), r2))

        L = len(aux1)
        for m in range(L):
            for n in range(m+1, L):
                if aux1[m][0] == aux1[n][0]:
                    angle = convfact * arccos(inner(aux1[m][3], aux1[n][3]) / (aux1[m][2] * aux1[n][2]))
                    ID.append(f'{aux1[m][1]} - {aux1[m][0]} - {aux1[n][1]}')
                    BA.append(angle)
                elif aux1[m][1] == aux1[n][1]:
                    angle = convfact * arccos(inner(aux1[m][3], aux1[n][3]) / (aux1[m][2] * aux1[n][2]))
                    ID.append(f'{aux1[m][0]} - {aux1[m][1]} - {aux1[n][0]}')
                    BA.append(angle)

    elif neigh1 == neigh2:
        nCenter = xyz.iloc[:,0].value_counts()[center]
        nNeigh = xyz.iloc[:,0].value_counts()[neigh1]

        xyzCenter = xyz[xyz['idAt'] == center]
        xyzNeigh = xyz[xyz['idAt'] == neigh1]

        rxCenter = xyzCenter.iloc[:, 1]
        ryCenter = xyzCenter.iloc[:, 2]
        rzCenter = xyzCenter.iloc[:, 3]

        rxNeigh = xyzNeigh.iloc[:, 1]
        ryNeigh = xyzNeigh.iloc[:, 2]
        rzNeigh = xyzNeigh.iloc[:, 3]

        for i in range(nCenter):
            for j in range(nNeigh):
                dx = rxNeigh.iloc[j] - rxCenter.iloc[i]
                dy = ryNeigh.iloc[j] - ryCenter.iloc[i]
                dz = rzNeigh.iloc[j] - rzCenter.iloc[i]

                dx = PBC(dx, Lx)
                dy = PBC(dy, Ly)
                dz = PBC(dz, Lz)
                r2 = [dx, dy, dz]

                d2 = dx**2 + dy**2 + dz**2
                if ((d2>= Rmin2) and (d2<= Rcut2)):
                    aux1.append((rxCenter.index[i],rxNeigh.index[j],sqrt(d2), r2))

        L = len(aux1)
        for m in range(L):
            for n in range(m+1, L):
                if aux1[m][0] == aux1[n][0]:
                    angle = convfact * arccos(inner(aux1[m][3], aux1[n][3]) / (aux1[m][2] * aux1[n][2]))
                    ID.append(f'{aux1[m][1]} - {aux1[m][0]} - {aux1[n][1]}')
                    BA.append(angle)

    elif neigh1 == center:
        nCenter = xyz.iloc[:,0].value_counts()[center]
        nNeigh = xyz.iloc[:,0].value_counts()[neigh2]

        xyzCenter = xyz[xyz['idAt'] == center]
        xyzNeigh = xyz[xyz['idAt'] == neigh2]

        rxCenter = xyzCenter.iloc[:, 1]
        ryCenter = xyzCenter.iloc[:, 2]
        rzCenter = xyzCenter.iloc[:, 3]

        rxNeigh = xyzNeigh.iloc[:, 1]
        ryNeigh = xyzNeigh.iloc[:, 2]
        rzNeigh = xyzNeigh.iloc[:, 3]

        for i in range(nCenter):
            for j in range(nNeigh):
                dx = rxNeigh.iloc[j] - rxCenter.iloc[i]
                dy = ryNeigh.iloc[j] - ryCenter.iloc[i]
                dz = rzNeigh.iloc[j] - rzCenter.iloc[i]

                dx = PBC(dx, Lx)
                dy = PBC(dy, Ly)
                dz = PBC(dz, Lz)
                r2 = [dx, dy, dz]

                d2 = dx**2 + dy**2 + dz**2
                if ((d2>= Rmin2) and (d2<= Rcut2)):
                    aux1.append((rxCenter.index[i],rxNeigh.index[j],sqrt(d2), r2))

        for i in range(nCenter):
            for j in range(i+1, nCenter):
                dx = rxCenter.iloc[j] - rxCenter.iloc[i]
                dy = ryCenter.iloc[j] - ryCenter.iloc[i]
                dz = rzCenter.iloc[j] - rzCenter.iloc[i]

                dx = PBC(dx, Lx)
                dy = PBC(dy, Ly)
                dz = PBC(dz, Lz)
                r2 = [dx, dy, dz]

                d2 = dx**2 + dy**2 + dz**2
                if ((d2>= Rmin2) and (d2<= Rcut2)):
                    aux2.append((rxCenter.index[i],rxCenter.index[j],sqrt(d2), r2))

        L1 = len(aux1)
        L2 = len(aux2)
        for m in range(L1):
            for n in range(L2):
                if aux1[m][0] == aux2[n][0]:
                    angle = convfact * arccos(inner(aux1[m][3], aux2[n][3]) / (aux1[m][2] * aux2[n][2]))
                    ID.append(f'{aux1[m][1]} - {aux1[m][0]} - {aux2[n][1]}')
                    BA.append(angle)
                elif aux1[m][0] == aux2[n][1]:
                    angle = convfact * arccos(inner(aux1[m][3], -aux2[n][3]) / (aux1[m][2] * aux2[n][2]))
                    ID.append(f'{aux1[m][1]} - {aux1[m][0]} - {aux2[n][0]}')
                    BA.append(angle)

    elif neigh2 == center:
        nCenter = xyz.iloc[:,0].value_counts()[center]
        nNeigh = xyz.iloc[:,0].value_counts()[neigh1]

        xyzCenter = xyz[xyz['idAt'] == center]
        xyzNeigh = xyz[xyz['idAt'] == neigh1]

        rxCenter = xyzCenter.iloc[:, 1]
        ryCenter = xyzCenter.iloc[:, 2]
        rzCenter = xyzCenter.iloc[:, 3]

        rxNeigh = xyzNeigh.iloc[:, 1]
        ryNeigh = xyzNeigh.iloc[:, 2]
        rzNeigh = xyzNeigh.iloc[:, 3]

        for i in range(nCenter):
            for j in range(nNeigh):
                dx = rxNeigh.iloc[j] - rxCenter.iloc[i]
                dy = ryNeigh.iloc[j] - ryCenter.iloc[i]
                dz = rzNeigh.iloc[j] - rzCenter.iloc[i]

                dx = PBC(dx, Lx)
                dy = PBC(dy, Ly)
                dz = PBC(dz, Lz)
                r2 = [dx, dy, dz]

                d2 = dx**2 + dy**2 + dz**2
                if ((d2>= Rmin2) and (d2<= Rcut2)):
                    aux1.append((rxCenter.index[i],rxNeigh.index[j],sqrt(d2), r2))

        for i in range(nCenter):
            for j in range(i+1, nCenter):
                dx = rxCenter.iloc[j] - rxCenter.iloc[i]
                dy = ryCenter.iloc[j] - ryCenter.iloc[i]
                dz = rzCenter.iloc[j] - rzCenter.iloc[i]

                dx = PBC(dx, Lx)
                dy = PBC(dy, Ly)
                dz = PBC(dz, Lz)
                r2 = [dx, dy, dz]

                d2 = dx**2 + dy**2 + dz**2
                if ((d2>= Rmin2) and (d2<= Rcut2)):
                    aux2.append((rxCenter.index[i],rxCenter.index[j],sqrt(d2), r2))

        L1 = len(aux1)
        L2 = len(aux2)
        for m in range(L1):
            for n in range(L2):
                if aux1[m][0] == aux2[n][0]:
                    angle = convfact * arccos(inner(aux1[m][3], aux2[n][3]) / (aux1[m][2] * aux2[n][2]))
                    ID.append(f'{aux1[m][1]} - {aux1[m][0]} - {aux2[n][1]}')
                    BA.append(angle)
                elif aux1[m][0] == aux2[n][1]:
                    angle = convfact * arccos(inner(aux1[m][3], -aux2[n][3]) / (aux1[m][2] * aux2[n][2]))
                    ID.append(f'{aux1[m][1]} - {aux1[m][0]} - {aux2[n][0]}')
                    BA.append(angle)

    else:
        nCenter = xyz.iloc[:,0].value_counts()[center]
        nNeigh1 = xyz.iloc[:,0].value_counts()[neigh1]
        nNeigh2 = xyz.iloc[:,0].value_counts()[neigh2]

        xyzCenter = xyz[xyz['idAt'] == center]
        xyzNeigh1 = xyz[xyz['idAt'] == neigh1]
        xyzNeigh2 = xyz[xyz['idAt'] == neigh2]

        rxCenter = xyzCenter.iloc[:, 1]
        ryCenter = xyzCenter.iloc[:, 2]
        rzCenter = xyzCenter.iloc[:, 3]

        rxNeigh1 = xyzNeigh1.iloc[:, 1]
        ryNeigh1 = xyzNeigh1.iloc[:, 2]
        rzNeigh1 = xyzNeigh1.iloc[:, 3]

        rxNeigh2 = xyzNeigh2.iloc[:, 1]
        ryNeigh2 = xyzNeigh2.iloc[:, 2]
        rzNeigh2 = xyzNeigh2.iloc[:, 3]

        for i in range(nCenter):
            for j in range(nNeigh1):
                dx = rxNeigh1.iloc[j] - rxCenter.iloc[i]
                dy = ryNeigh1.iloc[j] - ryCenter.iloc[i]
                dz = rzNeigh1.iloc[j] - rzCenter.iloc[i]

                dx = PBC(dx, Lx)
                dy = PBC(dy, Ly)
                dz = PBC(dz, Lz)
                r2 = [dx, dy, dz]

                d2 = dx**2 + dy**2 + dz**2
                if ((d2>= Rmin2) and (d2<= Rcut2)):
                    aux1.append((rxCenter.index[i],rxNeigh1.index[j],sqrt(d2), r2))

        for i in range(nCenter):
            for j in range(nNeigh2):
                dx = rxNeigh2.iloc[j] - rxCenter.iloc[i]
                dy = ryNeigh2.iloc[j] - ryCenter.iloc[i]
                dz = rzNeigh2.iloc[j] - rzCenter.iloc[i]

                dx = PBC(dx, Lx)
                dy = PBC(dy, Ly)
                dz = PBC(dz, Lz)
                r2 = [dx, dy, dz]

                d2 = dx**2 + dy**2 + dz**2
                if ((d2>= Rmin2) and (d2<= Rcut2)):
                    aux2.append((rxCenter.index[i],rxNeigh2.index[j],sqrt(d2), r2))

        L1 = len(aux1)
        L2 = len(aux2)
        for m in range(L1):
            for n in range(L2):
                if aux1[m][0] == aux2[n][0]:
                    angle = convfact * arccos(inner(aux1[m][3], aux2[n][3]) / (aux1[m][2] * aux2[n][2]))
                    ID.append(f'{aux1[m][1]} - {aux1[m][0]} - {aux2[n][1]}')
                    BA.append(angle)

    Summary = f'Summary >> {neigh1}-{center}-{neigh2} = ({mean(BA):.1f} +- {std(BA):.1f})° ; Count = {len(BA)} <<'

    with open(f'{neigh1}-{center}-{neigh2}.csv', 'w') as f:
        f.write('==== Bond Agnle in degrees ==== \n\n')
        f.write('Atoms ID, Angle \n')
        for i in range(len(BA)):
            f.write(f'{ID[i]}, {BA[i]:.2f} \n')
        f.write(f'\n {Summary}')

    print(f'Job done in {(time() - start):.3f} seconds!')
    print(f'Output file: {neigh1}-{center}-{neigh2}.csv')
    print(f'{Summary}')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('input_file', help = "Path to the xyz input file.")

    parser.add_argument('Rcut', type = float, help = "Maximum distance to be \
                        considered as bond length.")

    parser.add_argument('central_atom')

    parser.add_argument('neighboring_atom', nargs = 2)

    parser.add_argument('--Rmin', type = float,
                        help = "Minimum distance to be considered as bond length.")

    args = parser.parse_args()

    main()
