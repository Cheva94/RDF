#!/usr/local/bin/python3.9

'''
    Calculation: Radial Distribution Function (RDF).
    Description: todo en Angstrom!
    Written by: Ignacio J. Chevallier-Boutell.
    Dated: July, 2021.
'''

import pandas as pd
import numpy as np

def init(xsf = 'example.xsf'):
    '''
    Takes an .xsf file and acquire:
        * nAt: number of atoms.
        * rows: number of rows with relevant information per snapshot.
        * idAt: list of elements' names.
        * steps: total number of snapshots.
        * cell: cell dimensions along x (Lx), y (Ly) and z (Lz).
        * snaps: dataframe with idAt and positions for the nAt in the steps snapshots
    '''

    nAt = pd.read_csv(xsf, header = None, delim_whitespace = True, skiprows = 7, nrows = 1).to_numpy()[0][0]

    rows = nAt + 2 # ver si puedo eliminar esas 2 filas extras!

    idAt = pd.read_csv(xsf, header = None, names = ['idAt'], delim_whitespace = True, skiprows = 8, usecols = [0], nrows = nAt).drop_duplicates()['idAt'].values.tolist()

    steps = pd.read_csv(xsf, header = None, delim_whitespace = True, nrows = 1).to_numpy()[0][1]

    cell = pd.read_csv(xsf, header = None, delim_whitespace = True, skiprows = 3, nrows = 3).to_numpy()

    Lx, Ly, Lz = cell[0][0], cell[1][1], cell[2][2]

    snaps = pd.read_csv(xsf, header = None, delim_whitespace = True, skiprows = 6, usecols = [0,1,2,3], names=['idAt', 'rx', 'ry', 'rz'])

    return nAt, rows, idAt, steps, Lx, Ly, Lz, snaps

def min_im(coord, length):
    '''
    Minimum image convention for PBC.
    '''

    if coord < -0.5 * length:
        coord += length
    elif coord > 0.5 * length:
        coord -= length

    return coord

def hist_up(Rcut, Lx, Ly, Lz, xyz, nAt, dr, RDF):
    '''
    Updates histogram information with current snapshot.
    '''

    Rcut2 = Rcut * Rcut

    label = np.array(xyz[:,0])
    rx = np.array(xyz[:,1])
    ry = np.array(xyz[:,2])
    rz = np.array(xyz[:,3])

    for i in range(nAt-1):
        for j in range(i+1, nAt):
            dx = rx[i] - rx[j]
            dx = min_im(dx, Lx)

            dy = ry[i] - ry[j]
            dy = min_im(dy, Ly)

            dz = rz[i] - rz[j]
            dz = min_im(dz, Lz)

            r2 = dx * dx + dy * dy + dz * dz

            if r2 <= Rcut2:
                bin = int(np.sqrt(r2)/dr)
                RDF[bin] += 2 # no sé bien por qué suma 2 y no 1 ! - dice que 2 es la contribución de las partículas

def normalize(nBin, dr, RDF, nRDF, nAt, f):
    '''
    Normalize de function.
    '''

    rho=0.85 # densidad!

    for i in range(nBin):
        r = dr * (i + 0.5) # Distance
        vBin = ( (i+1) * (i+1) * (i+1) - i*i*i ) * dr*dr*dr # volume between bin i and bin i+1
        nId = (4/3) * np.pi * vBin * rho # number of ideal particles in vBin
        RDF[i] /= nRDF * nAt * nId

        f.write(f'{dr*i:.3} \t {RDF[i]:.4} \n')

def main():

    nAt, rows, idAt, steps, Lx, Ly, Lz, snaps = init()

    Rcut = 10 # maximum radius to be considered
    dr = 0.1 # bin size
    nBin = int(Rcut/dr) + 1 # number of bins
    RDF = np.zeros(nBin) # initialize array of zeros for RDF
    nRDF = 0 # number of snapshots used

    with open(f'rdf.dat', 'w') as f:
        for i in range(steps):
        # for i in range(1):
            nRDF += 1
            start = i*rows + 2
            end = rows*(i+1)
            xyz = snaps.iloc[start:end].to_numpy()
            hist_up(Rcut, Lx, Ly, Lz, xyz, nAt, dr, RDF)
        
        RDF = normalize(nBin, dr, RDF, nRDF, nAt, f)

if __name__ == "__main__":
    main()
