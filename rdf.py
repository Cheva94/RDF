#!/usr/local/bin/python3.9

'''
    Calculation: Radial Distribution Function (RDF).
    Description:
    Written by: Ignacio J. Chevallier-Boutell.
    Dated: July, 2021.
'''

import pandas as pd
import numpy as np

def init(xsf = 'example.xsf'):

    N = pd.read_csv(xsf, header = None, delim_whitespace = True, skiprows = 7, nrows = 1).to_numpy()[0][0]

    rows = N + 2

    ID = pd.read_csv(xsf, header = None, names = ['ID'], delim_whitespace = True, skiprows = 8, usecols = [0], nrows = N).drop_duplicates()['ID'].values.tolist()

    steps = pd.read_csv(xsf, header = None, delim_whitespace = True, nrows = 1).to_numpy()[0][1]

    cell = pd.read_csv(xsf, header = None, delim_whitespace = True, skiprows = 3, nrows = 3).to_numpy()

    snaps = pd.read_csv(xsf, header = None, delim_whitespace = True, skiprows = 6, usecols = [0,1,2,3], names=['ID', 'rx', 'ry', 'rz'])

    ID = pd.read_csv(xsf, header = None, names = ['ID'], delim_whitespace = True, skiprows = 8, usecols = [0], nrows = N).drop_duplicates()['ID'].values.tolist()

    return N, rows, ID, steps, cell, snaps

def min_im(coord, length):

    if coord < -0.5 * length:
        coord += length
    elif coord >= 0.5 * length:
        coord -= length

    return coord

def hist_up(N, cell, xyz, dr, Rcut, H):

    Rcut2 = Rcut * Rcut

    Lx, Ly, Lz = cell[0][0], cell[1][1], cell[2][2]

    label = np.array(xyz[:,0])
    rx = np.array(xyz[:,1])
    ry = np.array(xyz[:,2])
    rz = np.array(xyz[:,3])

    for i in range(N-1):
        for j in range(i+1, N):
            dx = rx[i] - rx[j]
            dx = min_im(dx, Lx)

            dy = ry[i] - ry[j]
            dy = min_im(dy, Ly)

            dz = rz[i] - rz[j]
            dz = min_im(dz, Lz)

            r2 = dx * dx + dy * dy + dz * dz

            if r2 < Rcut2:
                bin = int(np.sqrt(r2)/dr)
                H[bin] += 2 # no sé bien por qué suma 2 y no 1

def normalize(dr, rho):
    for i in range(nB):
        r = dr * (i + 0.5) # Distance
        vb = ( (i+1) * (i+1) * (i+1) - i*i*i ) * dr*dr*dr # volume between bin i and i+1
        nid = (4/3) * np.pi * vb * rho

def main():

    N, rows, ID, steps, cell, snaps = init()
    Rcut = 10 # Angstrom
    dr = 0.1 # Angstrom
    nB = int(Rcut/dr) + 1 # Total number of bins
    H = np.zeros(nB, dtype=int)

    rho=0.85
    ngr = 0

    # for i in range(steps):
    for i in range(1):
        ngr += 1
        start = i*rows + 2
        end = rows*(i+1)
        xyz = snaps.iloc[start:end].to_numpy()

    # faltan cosas acá




if __name__ == "__main__":
    main()
