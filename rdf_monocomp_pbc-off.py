#!python3.9

'''
    Calculation: Radial Distribution Function (RDF).
    Description: using Angstrom.
    Written by: Ignacio J. Chevallier-Boutell.
    Dated: July, 2021.
'''

from pandas import read_csv
from numpy import zeros, sqrt, array, pi

def user_input():
    '''
    User queries.
    '''

    # System
    # file_in = input("Enter your xsf file: ")
    file_in = '~/tesis/AIMD/AIMD_Agua-13x@silica03-2x2_10ps_WaterO2Pt.xsf'
    xsf = read_csv(file_in, header = None, delim_whitespace = True, names=['idAt', 'rx', 'ry', 'rz', 'fx', 'fy', 'fz'])

    total_frames = int(xsf.iloc[0,1])
    Lx, Ly, Lz = float(xsf.iloc[3,0]), xsf.iloc[4,1], xsf.iloc[5,2]
    nAtTot = int(xsf.iloc[7,0])
    idAt = xsf.iloc[8:(nAtTot+8),0].drop_duplicates().tolist()
    xyz_all = xsf.iloc[6:,0:4].reset_index(drop=True)

    # print('>>> Information from the input file:')
    # print(f'\t * There are {total_frames} frames in the file.')
    # print(f'\t * Cell dimensions (in Angstrom): (x, y, z) = ({Lx:.2f}, {Ly:.2f}, {Lz:.2f}).')
    # print(f'\t * There are {nAtTot} atoms whithin the cell.')
    # print(f'\t * Atomic species: {", ".join(idAt)}.')

    # Atoms to compare
    # at = input("Choose one element of the list of atomic species: ")
    at = 'Pt'

    nAt = xsf.iloc[8:(nAtTot+8),0].value_counts()[at]

    # # Number of frames
    # frame_start = 0
    # frame_end = -1

    # Histogram parameters
    dr = 0.1 # increment
    Rcut = 10.0 # maximum radius to be considered (max Value of the histogram)

    # Output file_out
    file_out = 'RDF_PBC-off_watO-watO'

    return total_frames, nAtTot, xyz_all, at, nAt, dr, Rcut, file_out

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
    Updates the existing histogram.

    It's considered that "data" is squared.
    '''

    binIdx = int(sqrt(data)/increment)
    H[binIdx] += 2 # contribution of i and j particles

def sample(xyz, dr, Rcut, RDF, nAt):
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

def normalize(nAt, dr, nBin, frames_count, RDF, file_out):
    '''
    Determine the normalize RDF.
    '''

    prefact = 4 * pi * dr**3

    with open(f'{file_out}.dat', 'w') as f:
        for binIdx in range(nBin):
            # r = [(binIdx+0.5)*dr for binIdx in range(nBin)] # distance to half bin
            volShell = prefact * (binIdx + 0.5)**2
            RDF[binIdx] /= frames_count * volShell
            f.write(f'{RDF[binIdx]} \n')

def main():

    total_frames, nAtTot, xyz_all, at, nAt, dr, Rcut, file_out = user_input()

    nBin, Rcut, RDF = hist_init(dr, Rcut)

    frames_count = 0
    rows = nAtTot + 2
    for frame in range(total_frames):
        xyz = xyz_all.iloc[(frame*rows + 2):((frame+1)*rows) , :]
        xyz = xyz[xyz['idAt'] == at].to_numpy()
        sample(xyz, dr, Rcut, RDF, nAt)
        frames_count += 1

    normalize(nAt, dr, nBin, frames_count, RDF, file_out)

    # print(f'Job done! The RDF file is: {file_out}.')

if __name__ == "__main__":
    main()
