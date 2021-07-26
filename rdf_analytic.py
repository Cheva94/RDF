#!python3.9

'''
    Calculation: Radial Distribution Function (RDF).
    Description: using Angstrom.
    Written by: Ignacio J. Chevallier-Boutell.
    Dated: July, 2021.
'''

import pandas as pd
from numpy import sqrt, zeros, array
from analytic_vol import volShell

def user_input():
    '''
    User queries.
    '''

    # System
    # file_in = input("Enter your xsf file: ")
    # file_in = 'example.xsf'
    file_in = 'largo.xsf'
    xsf = pd.read_csv(file_in, header = None, delim_whitespace = True, names=['idAt', 'rx', 'ry', 'rz', 'fx', 'fy', 'fz'])

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
    # at1 = input("Choose one element of the list of atomic species: ")
    # at2 = input("Choose another element of the list of atomic species: ")
    at1 = 'Si'
    at2 = 'Pt'

    nAt1 = xsf.iloc[8:(nAtTot+8),0].value_counts()[at1]
    nAt2 = xsf.iloc[8:(nAtTot+8),0].value_counts()[at2]

    # xyz1 = xyz_all[xyz_all['idAt'] == at1].to_numpy()
    # xyz2 = xyz_all[xyz_all['idAt'] == at2].to_numpy()

    # # Number of frames
    # frame_start = 0
    # frame_end = -1

    # Histogram parameters
    dr = 0.1 # increment
    Rcut = 0.5 * sqrt(Lx**2 + Ly**2 + Lz**2) # maximum radius to be considered (max Value of the histogram)

    # Output file_out
    file_out = f'{file_in.split(".")[0]}_{at1}-{at2}_rdf_analytic'

    return total_frames, Lx, Ly, Lz, nAtTot, xyz_all, at1, at2, nAt1, nAt2, dr, Rcut, file_out

def hist_init(increment, maximum):
    '''
    Initialize the histogram.
    '''

    nBin = int(maximum/increment) + 1 # number of bins
    maximum = nBin * increment # adjust maximum
    H = zeros(nBin) # initialize array of zeros for global histogram
    OccBin = zeros(nBin) # count bins with at least on object

    return nBin, maximum, H, OccBin

def hist_up(data, increment, H):
    '''
    Updates the existing histogram.

    It's considered that "data" is squared.
    '''

    binIdx = int(sqrt(data)/increment)
    H[binIdx] += 2 # contribution of i and j particles

def sample(Lx, Ly, Lz, xyz1, xyz2, dr, Rcut, global_RDF, nAt1, nAt2, nBin, OccShell):
    Rcut2 = Rcut * Rcut

    rx1 = array(xyz1[:,1])
    ry1 = array(xyz1[:,2])
    rz1 = array(xyz1[:,3])

    rx2 = array(xyz2[:,1])
    ry2 = array(xyz2[:,2])
    rz2 = array(xyz2[:,3])

    # Updates the histogram
    for i in range(nAt1):
        local_RDF = zeros(nBin)
        for j in range(nAt2):
            dx = rx1[i] - rx2[j]
            dy = ry1[i] - ry2[j]
            dz = rz1[i] - rz2[j]

            d2 = dx**2 + dy**2 + dz**2

            if d2 <= Rcut2:
                hist_up(d2, dr, local_RDF)

        localBox = [-rx1[i], Lx - rx1[i], -ry1[i], Ly - ry1[i], -rz1[i], Lz - rz1[i]]
        for binIdx in range(nBin):
            volBin = volShell(binIdx * dr, (binIdx + 1) * dr, localBox)
            if volBin > 0: # revisar esto!
                local_RDF[binIdx] /= volBin
                OccShell[binIdx] += 2 # decia +1 pero por como lo tengo escrito es +2

        global_RDF += local_RDF

def main():

    total_frames, Lx, Ly, Lz, nAtTot, xyz_all, at1, at2, nAt1, nAt2, dr, Rcut, file_out = user_input()

    nBin, Rcut, global_RDF, OccShell = hist_init(dr, Rcut)

    time_RDF = zeros(nBin)
    frames_count = 0
    rows = nAtTot + 2
    for frame in range(total_frames):
        xyz = xyz_all.iloc[(frame*rows + 2):((frame+1)*rows) , :]
        xyz1 = xyz[xyz['idAt'] == at1].to_numpy()
        xyz2 = xyz[xyz['idAt'] == at2].to_numpy()

        sample(Lx, Ly, Lz, xyz1, xyz2, dr, Rcut, global_RDF, nAt1, nAt2, nBin, OccShell)

        # rho = (nAt1 + nAt2) / (Lx * Ly * Lz) # densidad ideal (uniforme)
        global_RDF *= (Lx * Ly * Lz) / (nAt1 + nAt2)
        for binIdx in range(nBin): # list comprehension?
            if OccShell[binIdx] != 0:
                global_RDF[binIdx] /= OccShell[binIdx]

        time_RDF += global_RDF
        frames_count += 1 # aun no lo use

    time_RDF /= frames_count

    with open(f'{file_out}.dat', 'w') as f:
        for binIdx in range(nBin):
            f.write(f'{time_RDF[binIdx]} \n')

    # print(f'Job done! The RDF file is: {file_out}.')

if __name__ == "__main__":
    main()
