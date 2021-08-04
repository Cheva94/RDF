#!/usr/bin/env python3.9

'''
    Calculation: 3D Radial Distribution Function (RDF).
    Description: Determines de 3D RDF when given a xsf file. The comparison can
                be made between the same kind of atom (monocomponent) or between
                different species (multicomponent). Periodic boundary conditions
                can be turn on. Everything is in Angstrom.
    Written by: Ignacio J. Chevallier-Boutell.
    Dated: August, 2021.
'''

import argparse
from core import *
from time import time

def main():
    start = time() # starting wall time

    dr = args.dr
    Rcut = args.Rcut
    frames_count = 0
    frame_start = args.frames[0]
    if frame_start != 0:
        frame_start -= 1

    nBin, Rcut, RDF = hist_init(0, Rcut, dr)

    if args.periodic_boundary_conditions:
        if args.monocomponent:
            print(f'Running 3D RDF between {at} and {at} with PBC.')

            at = args.monocomponent

            total_frames, Lx, Ly, Lz, nAtTot, nAt, xyz_all = user_file_mono(args.input_file, at)

            Lmin = 0.5 * min(Lx, Ly, Lz)
            if Rcut > Lmin:
                print(f'Cannot choose Rcut greater than {Lmin:.2f}.')
                print('This will be the new Rcut value.')
                Rcut = Lmin

            rows = nAtTot + 2

            frame_end = args.frames[1]
            if frame_end == -1:
                frame_end = total_frames

            for frame in range(frame_start, frame_end):
                xyz = xyz_all.iloc[(frame*rows + 2):((frame+1)*rows) , :]
                xyz = xyz[xyz['idAt'] == at].to_numpy()
                mono_on_sample3d(Lx, Ly, Lz, xyz, dr, Rcut, RDF, nAt)
                frames_count += 1

            output_file = args.output_file
            if output_file == None:
                output_file = f'RDF3D_{at}-{at}_dr-{dr}_Rcut-{Rcut:.1f}_PBC'

            normalize_on_mono3d(Lx, Ly, Lz, nAt, dr, nBin, frames_count, RDF,
                                output_file)

            elapsed = time() - start # elapsed wall time
            print(f'Job done in {elapsed:.3f} seconds!')
            print(f'Output file: {output_file}.dat')

        elif args.multicomponents:
            print(f'Running 3D RDF between {at1} and {at2} with PBC.')

            at1 = args.multicomponents[0]
            at2 = args.multicomponents[1]

            total_frames, Lx, Ly, Lz, nAtTot, nAt1, nAt2, xyz_all = user_file_multi(args.input_file, at1, at2)

            Lmin = 0.5 * min(Lx, Ly, Lz)
            if Rcut > Lmin:
                print(f'Cannot choose Rcut greater than {Lmin:.2f}.')
                print('This will be the new Rcut value.')
                Rcut = Lmin

            rows = nAtTot + 2

            frame_end = args.frames[1]
            if frame_end == -1:
                frame_end = total_frames

            for frame in range(frame_start, frame_end):
                xyz = xyz_all.iloc[(frame*rows + 2):((frame+1)*rows) , :]
                xyz1 = xyz[xyz['idAt'] == at1].to_numpy()
                xyz2 = xyz[xyz['idAt'] == at2].to_numpy()
                multi_on_sample3d(Lx, Ly, Lz, xyz1, xyz2, dr, Rcut, RDF, nAt1, nAt2)
                frames_count += 1

            output_file = args.output_file
            if output_file == None:
                output_file = f'RDF3D_{at1}-{at2}_dr-{dr}_Rcut-{Rcut:.1f}_PBC'

            normalize_on_multi3d(Lx, Ly, Lz, nAt1, nAt2, dr, nBin, frames_count,
                                RDF, output_file)

            elapsed = time() - start # elapsed wall time
            print(f'Job done in {elapsed:.3f} seconds!')
            print(f'Output file: {output_file}.dat')

        else:
            print('Must choose mono or multi, and select elements to compare.')

    else:
        if args.monocomponent:
            print(f'Running 3D RDF between {at} and {at} without PBC.')

            at = args.monocomponent

            total_frames, Lx, Ly, Lz, nAtTot, nAt, xyz_all = user_file_mono(args.input_file, at)

            rows = nAtTot + 2

            frame_end = args.frames[1]
            if frame_end == -1:
                frame_end = total_frames

            for frame in range(frame_start, frame_end):
                xyz = xyz_all.iloc[(frame*rows + 2):((frame+1)*rows) , :]
                xyz = xyz[xyz['idAt'] == at].to_numpy()
                mono_off_sample3d(xyz, dr, Rcut, RDF, nAt)
                frames_count += 1

            output_file = args.output_file
            if output_file == None:
                output_file = f'RDF3D_{at}-{at}_dr-{dr}_Rcut-{Rcut:.1f}'

            normalize_off3d(dr, nBin, frames_count, RDF, output_file)

            elapsed = time() - start # elapsed wall time
            print(f'Job done in {elapsed:.3f} seconds!')
            print(f'Output file: {output_file}.dat')

        elif args.multicomponents:
            print(f'Running 3D RDF between {at1} and {at2} without PBC.')

            at1 = args.multicomponents[0]
            at2 = args.multicomponents[1]

            total_frames, Lx, Ly, Lz, nAtTot, nAt1, nAt2, xyz_all = user_file_multi(args.input_file, at1, at2)

            rows = nAtTot + 2

            frame_end = args.frames[1]
            if frame_end == -1:
                frame_end = total_frames

            for frame in range(frame_start, frame_end):
                xyz = xyz_all.iloc[(frame*rows + 2):((frame+1)*rows) , :]
                xyz1 = xyz[xyz['idAt'] == at1].to_numpy()
                xyz2 = xyz[xyz['idAt'] == at2].to_numpy()
                multi_off_sample3d(xyz1, xyz2, dr, Rcut, RDF, nAt1, nAt2)
                frames_count += 1

            output_file = args.output_file
            if output_file == None:
                output_file = f'RDF3D_{at1}-{at2}_dr-{dr}_Rcut-{Rcut:.1f}'

            normalize_off3d(dr, nBin, frames_count, RDF, output_file)

            elapsed = time() - start # elapsed wall time
            print(f'Job done in {elapsed:.3f} seconds!')
            print(f'Output file: {output_file}.dat')

        else:
            print('Must choose mono or multi, and select elements to compare.')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('-pbc', '--periodic_boundary_conditions',
                        action = 'store_true', help = "Set PBC on.")

    parser.add_argument('-mono', '--monocomponent', help = "Comparison between \
                        the same kind of atom. One argument needed.")

    parser.add_argument('-multi', '--multicomponents', nargs = 2,
                        help = "Comparison between different species. Two \
                        arguments needed.")

    parser.add_argument('input_file', help = "Path to the xsf input file.")

    parser.add_argument('dr', type = float, help = "Increment to be considered.")

    parser.add_argument('Rcut', type = float, help = "Maximum radius to be \
                        considered.")

    parser.add_argument('-o', '--output_file', help = "Path to the output file. \
                        If not given, the default name will be used.")

    parser.add_argument('-f', '--frames', type = int, nargs = 2, default = [0, -1],
                        help = "Choose starting and ending frames to compute.")

    args = parser.parse_args()

    main()
