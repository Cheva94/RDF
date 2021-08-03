#!/usr/bin/env python3.9

'''
    Calculation: 2D Radial Distribution Function (RDF).
    Description: Determines de 2D RDF when given a xsf file. The comparison can
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
    frame_start = int(args.frames[0])
    if frame_start != 0:
        frame_start -= 1

    dh = args.dh
    Hmin = args.Hcut[0]
    Hmax = args.Hcut[1]

    nBin, Rcut, RDF = hist_init(0, Rcut, dr)

    if args.periodic_boundary_conditions:
        if args.monocomponent:

            at = args.monocomponent

            total_frames, Lx, Ly, Lz, nAtTot, nAt, xyz_all = user_file_mono(args.input_file, at)

            Lmin = 0.5 * min(Lx, Ly)
            if Rcut > Lmin:
                print(f'Cannot choose Rcut greater than {Lmin:.2f}. This will be the new Rcut value.')
                Rcut = Lmin

            rows = nAtTot + 2

            frame_end = int(args.frames[1])
            if frame_end == -1:
                frame_end = total_frames

            for frame in range(frame_start, frame_end):
                xyz = xyz_all.iloc[(frame*rows + 2):((frame+1)*rows) , :]
                xyz = xyz[(xyz['idAt'] == at) & (h - 0.5 * dh <= xyz['rz']) & (xyz['rz'] <= h + 0.5 * dh)].to_numpy()
                nAt = len(xyz)
                mono_on_sample2d(Lx, Ly, xyz, dr, Rcut, RDF, nAt)
                frames_count += 1

            output_file = args.output_file
            if output_file == None:
                output_file = f'RDF3D_{at}-{at}_dr-{dr}_Rcut-{Rcut}_dh-{dh}_PBC'

            normalize_on_mono2d(Lx, Ly, dh, nAt, dr, nBin, frames_count, RDF,
                                output_file)

            elapsed = time() - start # elapsed wall time
            print(f'Job done in {elapsed:.3f} seconds!')
            print(f'2D RDF between {at} and {at} was calculated with PBC.')
            print(f'Output file: {output_file}.dat')

        elif args.multicomponents:

            at1 = args.multicomponents[0]
            at2 = args.multicomponents[1]

            total_frames, Lx, Ly, Lz, nAtTot, nAt1, nAt2, xyz_all = user_file_multi(args.input_file, at1, at2)

            Lmin = 0.5 * min(Lx, Ly)
            if Rcut > Lmin:
                print(f'Cannot choose Rcut greater than {Lmin:.2f}. This will be the new Rcut value.')
                Rcut = Lmin

            rows = nAtTot + 2

            frame_end = int(args.frames[1])
            if frame_end == -1:
                frame_end = total_frames

            for frame in range(frame_start, frame_end):
                xyz = xyz_all.iloc[(frame*rows + 2):((frame+1)*rows) , :]
                xyz1 = xyz[(xyz['idAt'] == at1) & (h - 0.5 * dh <= xyz['rz']) & (xyz['rz'] <= h + 0.5 * dh)].to_numpy()
                xyz2 = xyz[(xyz['idAt'] == at2) & (h - 0.5 * dh <= xyz['rz']) & (xyz['rz'] <= h + 0.5 * dh)].to_numpy()
                nAt1 = len(xyz1)
                nAt2 = len(xyz2)
                multi_on_sample2d(Lx, Ly, xyz1, xyz2, dr, Rcut, RDF, nAt1, nAt2)
                frames_count += 1

            output_file = args.output_file
            if output_file == None:
                output_file = f'RDF3D_{at1}-{at2}_dr-{dr}_Rcut-{Rcut}_dh-{dh}_PBC'

            normalize_on_multi2d(Lx, Ly, dh, nAt1, nAt2, dr, nBin, frames_count,
                                RDF, output_file)

            elapsed = time() - start # elapsed wall time
            print(f'Job done in {elapsed:.3f} seconds!')
            print(f'2D RDF between {at1} and {at2} was calculated with PBC.')
            print(f'Output file: {output_file}.dat')

        else:
            print('Must choose mono or multi, and select elements to compare.')

    else:
        if args.monocomponent:

            at = args.monocomponent

            total_frames, Lx, Ly, Lz, nAtTot, nAt, xyz_all = user_file_mono(args.input_file, at)

            rows = nAtTot + 2

            frame_end = int(args.frames[1])
            if frame_end == -1:
                frame_end = total_frames

            ##########################################
            nSlabs = int((Hmax - Hmin)/dh) + 1
            Hmax = nSlabs * dh + Hmin

            for slabIdx in range(nSlabs):
                h = (slabIdx + 0.5) * dh + Hmin

                for frame in range(frame_start, frame_end):
                    xyz = xyz_all.iloc[(frame*rows + 2):((frame+1)*rows) , :]
                    xyz = xyz[(xyz['idAt'] == at) & (h - 0.5 * dh <= xyz['rz']) & (xyz['rz'] <= h + 0.5 * dh)].to_numpy()
                    nAt = len(xyz)
                    mono_off_sample2d(xyz, dr, Rcut, RDF, nAt)
                    frames_count += 1

                output_file = args.output_file
                if output_file == None:
                    output_file = f'RDF3D_{at}-{at}_dr-{dr}_Rcut-{Rcut}_dh-{dh}'

                normalize_off2d(dh, dr, nBin, frames_count, RDF, output_file)

                elapsed = time() - start # elapsed wall time
                print(f'Job done in {elapsed:.3f} seconds!')
                print(f'2D RDF between {at} and {at} was calculated without PBC.')
                print(f'Output file: {output_file}.dat')

                RDF = zeros(nBin) # initialize array of zeros

            ##########################################

        elif args.multicomponents:


            at1 = args.multicomponents[0]
            at2 = args.multicomponents[1]

            total_frames, Lx, Ly, Lz, nAtTot, nAt1, nAt2, xyz_all = user_file_multi(args.input_file, at1, at2)

            rows = nAtTot + 2


            frame_end = int(args.frames[1])
            if frame_end == -1:
                frame_end = total_frames

            for frame in range(frame_start, frame_end):
                xyz = xyz_all.iloc[(frame*rows + 2):((frame+1)*rows) , :]
                xyz1 = xyz[(xyz['idAt'] == at1) & (h - 0.5 * dh <= xyz['rz']) & (xyz['rz'] <= h + 0.5 * dh)].to_numpy()
                xyz2 = xyz[(xyz['idAt'] == at2) & (h - 0.5 * dh <= xyz['rz']) & (xyz['rz'] <= h + 0.5 * dh)].to_numpy()
                nAt1 = len(xyz1)
                nAt2 = len(xyz2)
                multi_off_sample2d(xyz1, xyz2, dr, Rcut, RDF, nAt1, nAt2)
                frames_count += 1

            output_file = args.output_file
            if output_file == None:
                output_file = f'RDF3D_{at1}-{at2}_dr-{dr}_Rcut-{Rcut}_dh-{dh}'

            normalize_off2d(dh, dr, nBin, frames_count, RDF, output_file)

            elapsed = time() - start # elapsed wall time
            print(f'Job done in {elapsed:.3f} seconds!')
            print(f'2D RDF between {at1} and {at2} was calculated without PBC.')
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

    parser.add_argument('dr', type = float, help = "Increment to be considered \
                        for the RDF.")

    parser.add_argument('Rcut', type = float, help = "Maximum radius to be \
                        considered in the RDF.")

    parser.add_argument('dh', type = float, help = "Width of the slab to be \
                        considered")

    parser.add_argument('Hcut', type = float, nargs = 2, default = [0, -1],
                        help = "Minimum and maximum heights to be considered.")

    parser.add_argument('-o', '--output_file', help = "Path to the output file. \
                        If not given, the default name will be used.")

    parser.add_argument('-f', '--frames', nargs = 2, default = [0, -1],
                        help = "Choose starting and ending frames to compute.")

    args = parser.parse_args()

    main()
