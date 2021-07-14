#!python3.9

'''
    Calculation: Radial pair Distribution Function (RDF).
    Description: todo en Angstrom!
    Written by: Ignacio J. Chevallier-Boutell.
    Dated: July, 2021.
'''

import pandas as pd
import numpy as np

def user_input():
    '''
    User queries.
    '''
    # rows = nAt + 2 # ver si puedo eliminar esas 2 filas extras!
    #
    # snaps = pd.read_csv(xsf, header = None, delim_whitespace = True, skiprows = 6, usecols = [0,1,2,3], names=['idAt', 'rx', 'ry', 'rz'])
    # # snaps es solo un snap o sea solo un xyz acá


    # System
    # xsf = 'example1.xsf'
    #
    # total_frames = pd.read_csv(xsf, header = None, delim_whitespace = True, nrows = 1).to_numpy()[0][1]
    #
    # cell = pd.read_csv(xsf, header = None, delim_whitespace = True, skiprows = 3, nrows = 3).to_numpy()
    #
    # Lx, Ly, Lz = cell[0][0], cell[1][1], cell[2][2]
    #
    # nAtTot = pd.read_csv(xsf, header = None, delim_whitespace = True, skiprows = 7, nrows = 1).to_numpy()[0][0]
    #
    # idAt = pd.read_csv(xsf, header = None, names = ['idAt'], delim_whitespace = True, skiprows = 8, usecols = [0], nrows = nAt).drop_duplicates()['idAt'].values.tolist()


    ####################
    xyz = pd.read_csv('example1_xyz.xsf', header = None, delim_whitespace = True, skiprows = 2, usecols = [0,1,2,3], names=['idAt', 'rx', 'ry', 'rz'])
    # nAtTot = 103
    # idAt = ['O', 'Si', 'Pt', 'H']
    Lx, Ly, Lz = 10.5599958569, 14.9399970185, 21.9999878486

    # Atoms to compare
    at1 = 'Si'
    at2 = 'Pt'

    xyz1 = xyz[xyz['idAt'] == at1].to_numpy()
    xyz2 = xyz[xyz['idAt'] == at2].to_numpy()

    # Number of frames

    # Histogram parameters
    dr = 0.1 # increment
    Rcut = 10.0 # maximum radius to be considered (max Value of the histogram)

    # Output filename
    filename = 'example1_NC'

    return Lx, Ly, Lz, xyz1, xyz2, dr, Rcut, filename

def hist_init(dr, Rcut):
    '''
    Initialize the histogram.
    '''

    nBin = int(Rcut/dr) + 1 # number of bins
    Rcut = nBin * dr # adjust maximum
    H = np.zeros(nBin) # initialize array of zeros

    return nBin, Rcut, H

def PBC(dist, length): # Nacho
    '''
    Correct a distance using Periodic Boundary Conditions (PBC).
    '''

    if dist < -0.5 * length:
        dist += length
    elif dist > 0.5 * length:
        dist -= length

    return dist

# def PBC(dist, length): # Frenkel
#     '''
#     Correct a distance using Periodic Boundary Conditions (PBC).
#     '''
#
#     return dist - length * int(dist/length)

def hist_up(data, dr, H):
    '''
    Updates the existing histogram.

    It's considered that "data" is squared.
    '''

    binIdx = int(np.sqrt(data)/dr)
    H[binIdx] += 2 # contribution of i and j particles

def sample(Lx, Ly, Lz, xyz1, xyz2, dr, Rcut, RDF):
    Rcut2 = Rcut * Rcut

    rx1 = np.array(xyz1[:,1])
    ry1 = np.array(xyz1[:,2])
    rz1 = np.array(xyz1[:,3])

    rx2 = np.array(xyz2[:,1])
    ry2 = np.array(xyz2[:,2])
    rz2 = np.array(xyz2[:,3])

    nAt1 = len(rx1)
    nAt2 = len(rx2)

    # Apply PBC and updates the histogram
    # No estoy considerando duplicados si midiera el mismo átomo contra sí mismo
    # Frenkel usa las PBC de manera diferente
    for i in range(nAt1):
        for j in range(nAt2):
            dx = rx1[i] - rx2[j]
            dx = PBC(dx, Lx)

            dy = ry1[i] - ry2[j]
            dy = PBC(dy, Ly)

            dz = rz1[i] - rz2[j]
            dz = PBC(dz, Lz)

            r2 = dx * dx + dy * dy + dz * dz

            if r2 <= Rcut2:
                # acá tiene 0.5 L siendo L el largo de la caja # Ver si va a ser contado y contar
                hist_up(r2, dr, RDF)

    return nAt1, nAt2

def normalize(Lx, Ly, Lz, nAt1, nAt2, dr, nBin, frames, RDF, filename):
    '''
    Determine the normalize RDF.
    '''

    nAt = nAt1 + nAt2 # number of atoms compared
    rho = nAt / (Lx * Ly * Lz) # ideal (uniform) density
    prefact = 4 * np.pi
    dr3 = dr * dr * dr

    with open(f'{filename}.dat', 'w') as f:
        for binIdx in range(nBin):
            # the distance is r = dr * (i + 0.5)
            volBin = prefact * ( (binIdx + 0.5) * (binIdx + 0.5) ) * dr3
            nIdeal = volBin * rho
            RDF[binIdx] /= frames * nIdeal * nAt
            f.write(f'{RDF[binIdx]} \n')

def main():

    Lx, Ly, Lz, xyz1, xyz2, dr, Rcut, filename = user_input()

    nBin, Rcut, RDF = hist_init(dr, Rcut)

    frames = 0 # frames counter

    nAt1, nAt2 = sample(Lx, Ly, Lz, xyz1, xyz2, dr, Rcut, RDF)
    frames += 1

    normalize(Lx, Ly, Lz, nAt1, nAt2, dr, nBin, frames, RDF, filename)

if __name__ == "__main__":
    main()
