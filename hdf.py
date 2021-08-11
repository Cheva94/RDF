 #!/usr/bin/env python3.9

'''
    Calculation: Height Distribution Function (HDF).
    Description: Determines the HDF along Z axis when given a xsf file.
                Everything is in Angstrom.
    Written by: Ignacio J. Chevallier-Boutell.
    Dated: August, 2021.
'''

import argparse
from coreHP import *
from time import time

def main():
    start = time()

    dh = args.dh
    Hmin = args.Hcut[0]
    Hmax = args.Hcut[1]
    at = args.at

    frames_count = 0
    frames_start = args.frames[0]
    if frames_start != 0:
        frames_start -= 1

    print(f'Running HDF for {at} atoms.')

    frames_total, Lx, Ly, Lz, nAtTot, nAt, xyz_all = userfile_mono(args.input_file, at)

    if Hmax == -1:
        Hmax = Lz

    nBin, Hmax, HDF = hist_init_hdf(Hmin, Hmax, dh)

    rows = nAtTot + 2

    frames_end = args.frames[1]
    if frames_end == -1:
        frames_end = frames_total

    for frame in range(frames_start, frames_end):
        xyz = xyz_all.iloc[(frame * rows + 2) : ((frame + 1) * rows), :]
        xyz = xyz[(xyz['idAt'] == at) & (Hmin <= xyz['rz']) & (xyz['rz'] < Hmax)].to_numpy()
        nAt = len(xyz)
        sample_hdf(xyz, dh, HDF, nAt, Hmin)
        frames_count += 1

    output_file = args.output_file
    if output_file == None:
        output_file = f'HDF:{at}_H:{Hmin}-{Hmax:.1f}_dh:{dh}'

    normalize_hdf(dh, nBin, frames_count, HDF, output_file, Hmin)

    print(f'Job done in {(time() - start):.3f} seconds!')
    print(f'Output file: {output_file}.csv')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('input_file', help = "Path to the xsf input file.")

    parser.add_argument('at', help = "Atom to be analyzed.")

    parser.add_argument('dh', type = float, help = "Increment to be considered.")

    parser.add_argument('-f', '--frames', type = int, nargs = 2, default = [0, -1],
                        help = "Choose starting and ending frames to compute.")

    parser.add_argument('-o', '--output_file', help = "Path to the output file. \
                        If not given, the default name will be used.")

    parser.add_argument('-H', '--Hcut', type = float, nargs = 2, default = [0, -1],
                        help = "Minimum and maximum heights to be considered.")

    args = parser.parse_args()

    main()
