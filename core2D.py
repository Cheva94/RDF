#!/usr/bin/env python3.9

'''
    Description: Core functions used in 2D RDF.
    Written by: Ignacio J. Chevallier-Boutell.
    Dated: August, 2021.
'''

from pandas import read_csv
from numpy import zeros, sqrt, array, pi

################################################################################
######################## Input processing functions
################################################################################

def userfile_mono(input_file, atom):
    '''
    Process the xsf input file given by the user in the monocomponent case.
    '''

    xsf = read_csv(input_file, header = None, delim_whitespace = True,
                    names=['idAt', 'rx', 'ry', 'rz', 'fx', 'fy', 'fz'])

    frames_total = int(xsf.iloc[0,1])
    Lx, Ly, Lz = float(xsf.iloc[3,0]), xsf.iloc[4,1], xsf.iloc[5,2]
    nAtTot = int(xsf.iloc[7,0])
    nAt = xsf.iloc[8:(nAtTot+8),0].value_counts()[atom]
    xyz_all = xsf.iloc[6:,0:4].reset_index(drop=True)

    return frames_total, Lx, Ly, Lz, nAtTot, nAt, xyz_all

def userfile_multi(input_file, atom1, atom2):
    '''
    Process the xsf input file given by the user in the multicomponent case.
    '''

    xsf = read_csv(input_file, header = None, delim_whitespace = True,
                    names=['idAt', 'rx', 'ry', 'rz', 'fx', 'fy', 'fz'])

    frames_total = int(xsf.iloc[0,1])
    Lx, Ly, Lz = float(xsf.iloc[3,0]), xsf.iloc[4,1], xsf.iloc[5,2]
    nAtTot = int(xsf.iloc[7,0])
    nAt1 = xsf.iloc[8:(nAtTot+8),0].value_counts()[atom1]
    nAt2 = xsf.iloc[8:(nAtTot+8),0].value_counts()[atom2]
    xyz_all = xsf.iloc[6:,0:4].reset_index(drop=True)

    return frames_total, Lx, Ly, Lz, nAtTot, nAt1, nAt2, xyz_all

################################################################################
######################## Histogramming functions
################################################################################

def hist_init(minimum, maximum, increment):
    '''
    Initialize 2D histogram.
    '''

    nBin = int((maximum - minimum)/increment) + 1
    maximum = nBin * increment + minimum
    H = zeros(nBin)

    return nBin, maximum, H

def hist_up(data, increment, H):
    '''
    Updates the existing 2D histogram.
    '''

    binIdx = int(data/increment)
    H[binIdx] += 1

################################################################################
######################## Sampling functions
################################################################################

def sample_off_mono(xyz, dr, Rcut, RDF, nAt):
    '''
    Determines 2D RDF in the monocomponent case without PBC.
    '''

    Rcut2 = Rcut * Rcut

    rx1 = array(xyz[:,1])
    ry1 = array(xyz[:,2])

    rx2 = array(xyz[:,1])
    ry2 = array(xyz[:,2])

    for i in range(nAt):
        for j in range(i+1, nAt):
            dx = rx1[i] - rx2[j]
            dy = ry1[i] - ry2[j]

            d2 = dx**2 + dy**2

            if d2 <= Rcut2:
                hist_up(sqrt(d2), dr, RDF)

def sample_off_multi(xyz1, xyz2, dr, Rcut, RDF, nAt1, nAt2):
    '''
    Determines 2D RDF in the multicomponent case without PBC.
    '''

    Rcut2 = Rcut * Rcut

    rx1 = array(xyz1[:,1])
    ry1 = array(xyz1[:,2])

    rx2 = array(xyz2[:,1])
    ry2 = array(xyz2[:,2])

    for i in range(nAt1):
        for j in range(nAt2):
            dx = rx1[i] - rx2[j]
            dy = ry1[i] - ry2[j]

            d2 = dx**2 + dy**2

            if d2 <= Rcut2:
                hist_up(sqrt(d2), dr, RDF)

def PBC(dist, length):
    '''
    Correct a distance using Periodic Boundary Conditions (PBC).
    '''

    return dist - length * int(2*dist/length)

def sample_on_mono(Lx, Ly, xyz, dr, Rcut, RDF, nAt):
    '''
    Determines 2D RDF in the monocomponent case with PBC.
    '''

    Rcut2 = Rcut * Rcut

    rx1 = array(xyz[:,1])
    ry1 = array(xyz[:,2])

    rx2 = array(xyz[:,1])
    ry2 = array(xyz[:,2])

    for i in range(nAt):
        for j in range(i+1, nAt):
            dx = rx1[i] - rx2[j]
            dy = ry1[i] - ry2[j]

            dx = PBC(dx, Lx)
            dy = PBC(dy, Ly)

            d2 = dx**2 + dy**2

            if d2 <= Rcut2:
                hist_up(sqrt(d2), dr, RDF)

def sample_on_multi(Lx, Ly, xyz1, xyz2, dr, Rcut, RDF, nAt1, nAt2):
    '''
    Determines 2D RDF in the multicomponent case with PBC.
    '''

    Rcut2 = Rcut * Rcut

    rx1 = array(xyz1[:,1])
    ry1 = array(xyz1[:,2])

    rx2 = array(xyz2[:,1])
    ry2 = array(xyz2[:,2])

    for i in range(nAt1):
        for j in range(nAt2):
            dx = rx1[i] - rx2[j]
            dy = ry1[i] - ry2[j]

            dx = PBC(dx, Lx)
            dy = PBC(dy, Ly)

            d2 = dx**2 + dy**2

            if d2 <= Rcut2:
                hist_up(sqrt(d2), dr, RDF)

################################################################################
######################## Normalizing functions
################################################################################

def normalize_off(dr, nBin, frames_count, RDF, output_file):
    '''
    Normalize the 2D RDF without PBC.
    '''

    prefact = pi * dr**2

    RDF /= 0.5 * frames_count

    with open(f'{output_file}.csv', 'w') as f:
        for binIdx in range(nBin):
            r = (binIdx + 0.5) * dr
            surfRing = prefact * (2 * binIdx + 1)
            RDF[binIdx] /= surfRing
            f.write(f'{r:.2f}, {RDF[binIdx]:.4f} \n')

def normalize_on_mono(Lx, Ly, dh, nAt, dr, nBin, frames_count, RDF, output_file):
    '''
    Normalize the 2D RDF with PBC in the monocomponent case.
    '''

    prefact = pi * dr**2
    volBox = Lx * Ly * dh
    nPairs = nAt * (nAt - 1)

    RDF *= (2 * volBox) / (nPairs * frames_count)

    with open(f'{output_file}.csv', 'w') as f:
        for binIdx in range(nBin):
            r = (binIdx + 0.5) * dr
            surfRing = prefact * (2 * binIdx + 1)
            RDF[binIdx] /= surfRing
            f.write(f'{r:.2f}, {RDF[binIdx]:.4f} \n')

def normalize_on_multi(Lx, Ly, dh, nAt1, nAt2, dr, nBin, frames_count, RDF,
                        output_file):
    '''
    Normalize the 2D RDF with PBC in the multicomponent case.
    '''

    prefact = pi * dr**2
    volBox = Lx * Ly * dh
    nPairs = nAt1 * nAt2 * 2

    RDF *= (2 * volBox) / (nPairs * frames_count)

    with open(f'{output_file}.csv', 'w') as f:
        for binIdx in range(nBin):
            r = (binIdx + 0.5) * dr
            surfRing = prefact * (2 * binIdx + 1)
            RDF[binIdx] /= surfRing
            f.write(f'{r:.2f}, {RDF[binIdx]:.4f} \n')
