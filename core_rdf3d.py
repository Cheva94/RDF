#!/usr/bin/env python3.9

'''
    Description: Core functions used in program rdf3d.py.
    Written by: Ignacio J. Chevallier-Boutell.
    Dated: July, 2021.
'''

from pandas import read_csv
from numpy import zeros, sqrt, array, pi

def user_file_mono(input_file, atom):
    '''
    Process the input file given by the user with -mono option.
    '''

    xsf = read_csv(input_file, header = None, delim_whitespace = True,
                    names=['idAt', 'rx', 'ry', 'rz', 'fx', 'fy', 'fz'])

    total_frames = int(xsf.iloc[0,1])
    Lx, Ly, Lz = float(xsf.iloc[3,0]), xsf.iloc[4,1], xsf.iloc[5,2]
    nAtTot = int(xsf.iloc[7,0])
    nAt = xsf.iloc[8:(nAtTot+8),0].value_counts()[atom]
    xyz_all = xsf.iloc[6:,0:4].reset_index(drop=True)

    return total_frames, Lx, Ly, Lz, nAtTot, nAt, xyz_all

def user_file_multi(input_file, atom1, atom2):
    '''
    Process the input file given by the user with -multi option.
    '''

    xsf = read_csv(input_file, header = None, delim_whitespace = True,
                    names=['idAt', 'rx', 'ry', 'rz', 'fx', 'fy', 'fz'])

    total_frames = int(xsf.iloc[0,1])
    Lx, Ly, Lz = float(xsf.iloc[3,0]), xsf.iloc[4,1], xsf.iloc[5,2]
    nAtTot = int(xsf.iloc[7,0])
    nAt1 = xsf.iloc[8:(nAtTot+8),0].value_counts()[atom1]
    nAt2 = xsf.iloc[8:(nAtTot+8),0].value_counts()[atom2]
    xyz_all = xsf.iloc[6:,0:4].reset_index(drop=True)

    return total_frames, Lx, Ly, Lz, nAtTot, nAt1, nAt2, xyz_all

def hist_init(increment, maximum):
    '''
    Initialize the histogram.
    '''

    nBin = int(maximum/increment) + 1 # number of bins
    maximum = nBin * increment # adjust maximum
    H = zeros(nBin) # initialize array of zeros

    return nBin, maximum, H

def hist_up(data, increment, H):
    '''
    Updates the existing histogram. It's considered that "data" is squared.
    '''

    binIdx = int(sqrt(data)/increment)
    H[binIdx] += 2 # contribution of i and j particles

def PBC(dist, length):
    '''
    Correct a distance using Periodic Boundary Conditions (PBC).
    '''
    return dist - length * int(2*dist/length)

def multi_off_sample(xyz1, xyz2, dr, Rcut, RDF, nAt1, nAt2):
    '''
    Determines RDF with -multi option and PBC off.
    '''

    Rcut2 = Rcut * Rcut

    rx1 = array(xyz1[:,1])
    ry1 = array(xyz1[:,2])
    rz1 = array(xyz1[:,3])

    rx2 = array(xyz2[:,1])
    ry2 = array(xyz2[:,2])
    rz2 = array(xyz2[:,3])

    # Updates the histogram
    for i in range(nAt1):
        for j in range(nAt2):
            dx = rx1[i] - rx2[j]
            dy = ry1[i] - ry2[j]
            dz = rz1[i] - rz2[j]

            d2 = dx**2 + dy**2 + dz**2

            if d2 <= Rcut2:
                hist_up(d2, dr, RDF)

def multi_on_sample(Lx, Ly, Lz, xyz1, xyz2, dr, Rcut, RDF, nAt1, nAt2):
    '''
    Determines RDF with -multi option and PBC on.
    '''

    Rcut2 = Rcut * Rcut

    rx1 = array(xyz1[:,1])
    ry1 = array(xyz1[:,2])
    rz1 = array(xyz1[:,3])

    rx2 = array(xyz2[:,1])
    ry2 = array(xyz2[:,2])
    rz2 = array(xyz2[:,3])

    # Apply PBC and updates the histogram
    for i in range(nAt1):
        for j in range(nAt2):
            dx = rx1[i] - rx2[j]
            dy = ry1[i] - ry2[j]
            dz = rz1[i] - rz2[j]

            dx = PBC(dx, Lx)
            dy = PBC(dy, Ly)
            dz = PBC(dz, Lz)

            d2 = dx**2 + dy**2 + dz**2

            if d2 <= Rcut2:
                hist_up(d2, dr, RDF)

def mono_on_sample(Lx, Ly, Lz, xyz, dr, Rcut, RDF, nAt):
    '''
    Determines RDF with -mono option and PBC on.
    '''

    Rcut2 = Rcut * Rcut

    rx1 = array(xyz[:,1])
    ry1 = array(xyz[:,2])
    rz1 = array(xyz[:,3])

    rx2 = array(xyz[:,1])
    ry2 = array(xyz[:,2])
    rz2 = array(xyz[:,3])

    # Apply PBC and updates the histogram
    for i in range(nAt):
        for j in range(i+1, nAt):
            dx = rx1[i] - rx2[j]
            dy = ry1[i] - ry2[j]
            dz = rz1[i] - rz2[j]

            dx = PBC(dx, Lx)
            dy = PBC(dy, Ly)
            dz = PBC(dz, Lz)

            d2 = dx**2 + dy**2 + dz**2

            if d2 <= Rcut2:
                hist_up(d2, dr, RDF)

def mono_off_sample(xyz, dr, Rcut, RDF, nAt):
    '''
    Determines RDF with -mono option and PBC off.
    '''

    Rcut2 = Rcut * Rcut

    rx1 = array(xyz[:,1])
    ry1 = array(xyz[:,2])
    rz1 = array(xyz[:,3])

    rx2 = array(xyz[:,1])
    ry2 = array(xyz[:,2])
    rz2 = array(xyz[:,3])

    # Updates the histogram
    for i in range(nAt):
        for j in range(i+1, nAt):
            dx = rx1[i] - rx2[j]
            dy = ry1[i] - ry2[j]
            dz = rz1[i] - rz2[j]

            d2 = dx**2 + dy**2 + dz**2

            if d2 <= Rcut2:
                hist_up(d2, dr, RDF)

def mono_on_normalize(Lx, Ly, Lz, nAt, dr, nBin, frames_count, RDF, output_file):
    '''
    Normalize de RDF with -mono option and PBC on.
    '''

    volBox = Lx * Ly * Lz
    nPairs = nAt * (nAt - 1)
    RDF *= volBox / nPairs
    prefact = 4 * pi * dr**3

    with open(f'{output_file}.dat', 'w') as f:
        for binIdx in range(nBin):
            r = (binIdx + 0.5) * dr
            volShell = prefact * (binIdx + 0.5)**2
            RDF[binIdx] /= frames_count * volShell
            f.write(f'{r:.2f}, {RDF[binIdx]:.4f} \n')

def mono_off_normalize(nAt, dr, nBin, frames_count, RDF, output_file):
    '''
    Normalize de RDF with -mono option and PBC off.
    '''

    prefact = 4 * pi * dr**3

    with open(f'{output_file}.dat', 'w') as f:
        for binIdx in range(nBin):
            r = (binIdx + 0.5) * dr
            volShell = prefact * (binIdx + 0.5)**2
            RDF[binIdx] /= frames_count * volShell
            f.write(f'{r:.2f}, {RDF[binIdx]:.4f} \n')

def multi_on_normalize(Lx, Ly, Lz, nAt1, nAt2, dr, nBin, frames_count, RDF,
                        output_file):
    '''
    Normalize de RDF with -multi option and PBC on.
    '''

    volBox = Lx * Ly * Lz
    nPairs = nAt1 * nAt2 * 2
    RDF *= volBox / nPairs
    prefact = 4 * pi * dr**3

    with open(f'{output_file}.dat', 'w') as f:
        for binIdx in range(nBin):
            r = (binIdx + 0.5) * dr
            volShell = prefact * (binIdx + 0.5)**2
            RDF[binIdx] /= frames_count * volShell
            f.write(f'{r:.2f}, {RDF[binIdx]:.4f} \n')

def multi_off_normalize(nAt1, nAt2, dr, nBin, frames_count, RDF, output_file):
    '''
    Normalize de RDF with -multi option and PBC off.
    '''

    prefact = 4 * pi * dr**3

    with open(f'{output_file}.dat', 'w') as f:
        for binIdx in range(nBin):
            r = (binIdx + 0.5) * dr
            volShell = prefact * (binIdx + 0.5)**2
            RDF[binIdx] /= frames_count * volShell
            f.write(f'{r:.2f}, {RDF[binIdx]:.4f} \n')
