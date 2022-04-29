#!/usr/bin/python3.9

'''
    Calculation: Height Distribution Function (HDF).
    Description: Determines the HDF along Z axis when given a xsf file.
                 Everything is in Angstrom.
    Written by: Ignacio J. Chevallier-Boutell.
    Dated: , 2021.
    Last update: April, 2022.
'''

import argparse
from core.coreHP import *
from time import time, strftime

def main():
    start = time()
    started = strftime("%H:%M:%S")

    dh = args.dh
    Hmin = args.Hcut[0]
    Hmax = args.Hcut[1]
    atL = args.at

    frames_start = args.frames[0]
    if frames_start != 0:
        frames_start -= 1

    for at in atL:

        print('############################')
        print(f'Running HDF for {at} atoms. Started: {started}.')

        frames_total, Lx, Ly, Lz, nAtTot, nAt, xyz_all = userfile_mono(args.input_file, at)

        if Hmax == -1:
            Hmax = Lz

        nBin, histMax, HDF = hist_init_hdf(Hmin, Hmax, dh)

        rows = nAtTot + 2

        frames_end = args.frames[1]
        if frames_end == -1:
            frames_end = frames_total

        nAtSlab = 0
        frames_count = 0

        for frame in range(frames_start, frames_end):
            xyz = xyz_all.iloc[(frame * rows + 2) : ((frame + 1) * rows), :]
            xyz = xyz[(xyz['idAt'] == at) & (Hmin <= xyz['rz']) & (xyz['rz'] < Hmax)].to_numpy()
            nAt = len(xyz)
            sample_hdf(xyz, dh, HDF, nAt, Hmin)
            frames_count += 1
            nAtSlab += nAt

        nAtAve = nAtSlab/frames_count

        output_file = args.output_file
        if output_file == None:
            output_file = f'HDF_{at}x{nAtAve:.0f}'

        normalize_hdf(dh, nBin, frames_count, HDF, output_file, Hmin)

        print(f'Job done in {(time() - start)/60:.3f} minutes!')
        print(f'Output file: {output_file}.csv')
        print(f'There are {nAtAve:.2f} {at} atoms on average within this slab.')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('input_file', help = "Path to the xsf input file.")

    parser.add_argument('dh', type = float, help = "Increment to be considered.")

    parser.add_argument('at', help = "Atom to be analyzed.", nargs = "+")

    parser.add_argument('-f', '--frames', type = int, nargs = 2, default = [0, -1],
                        help = "Choose starting and ending frames to compute.")

    parser.add_argument('-o', '--output_file', help = "Path to the output file. \
                        If not given, the default name will be used.")

    parser.add_argument('-H', '--Hcut', type = float, nargs = 2, default = [0, -1],
                        help = "Minimum and maximum heights to be considered.")

    args = parser.parse_args()

    main()
