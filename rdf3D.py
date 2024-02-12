#!/usr/bin/python3.10

'''
    Calculation: 3D Radial Distribution Function (RDF).
    Description: Determines the 3D RDF when given a xsf file. The comparison can
                be made between the same kind of atom (monocomponent) or between
                different species (multicomponent). Periodic boundary conditions
                can be turn on with -pbc option. Everything is in Angstrom.
    Written by: Ignacio J. Chevallier-Boutell.
    Dated: August, 2021.
'''

import argparse
from core.core3D import *
from time import time, strftime

def main():
    start = time()
    started = strftime("%H:%M:%S")

    Rcut = args.Rcut
    dr = args.dr
    frames_count = 0
    frames_start = args.frames[0]
    if frames_start != 0:
        frames_start -= 1

    if args.periodic_boundary_conditions:
        if args.monocomponent:
            at = args.monocomponent
            print(f'Running 3D RDF between {at} atoms with PBC. Started: {started}.')

            frames_total, Lx, Ly, Lz, nAtTot, nAt, xyz_all = userfile_mono(args.input_file, at)

            Lmin = 0.5 * min(Lx, Ly, Lz)
            if Rcut > Lmin:
                print(f'Cannot choose Rcut greater than {Lmin:.3f}.')
                print('This will be the new Rcut value.')
                Rcut = Lmin

            nBin, Rcut, RDF = hist_init(Rcut, dr)

            rows = nAtTot + 2

            frames_end = args.frames[1]
            if frames_end == -1:
                frames_end = frames_total

            for frame in range(frames_start, frames_end):
                xyz = xyz_all.iloc[(frame * rows + 2) : ((frame + 1) * rows), :]
                xyz = xyz[xyz['idAt'] == at].to_numpy()
                sample_on_mono(Lx, Ly, Lz, xyz, dr, Rcut, RDF, nAt)
                frames_count += 1

            output_file = args.output_file
            if output_file == None:
                output_file = f'RDF3D_{at}-{at}_PBC'

            normalize_on_mono(Lx, Ly, Lz, nAt, dr, nBin, frames_count, RDF,
                                output_file)

            print(f'Job done in {(time() - start)/60:.3f} minutes!')
            print(f'Output file: {output_file}.csv')

        elif args.multicomponents:
            at1 = args.multicomponents[0]
            at2 = args.multicomponents[1]
            print(f'Running 3D RDF between {at1} and {at2} atoms with PBC. Started: {started}.')

            frames_total, Lx, Ly, Lz, nAtTot, nAt1, nAt2, xyz_all = userfile_multi(args.input_file, at1, at2)

            Lmin = 0.5 * min(Lx, Ly, Lz)
            if Rcut > Lmin:
                print(f'Cannot choose Rcut greater than {Lmin:.3f}.')
                print('This will be the new Rcut value.')
                Rcut = Lmin

            nBin, Rcut, RDF = hist_init(Rcut, dr)

            rows = nAtTot + 2

            frames_end = args.frames[1]
            if frames_end == -1:
                frames_end = frames_total

            for frame in range(frames_start, frames_end):
                xyz = xyz_all.iloc[(frame * rows + 2) : ((frame + 1) * rows), :]
                xyz1 = xyz[xyz['idAt'] == at1].to_numpy()
                xyz2 = xyz[xyz['idAt'] == at2].to_numpy()
                sample_on_multi(Lx, Ly, Lz, xyz1, xyz2, dr, Rcut, RDF, nAt1, nAt2)
                frames_count += 1

            output_file = args.output_file
            if output_file == None:
                output_file = f'RDF3D_{at1}-{at2}_PBC'

            normalize_on_multi(Lx, Ly, Lz, nAt1, nAt2, dr, nBin, frames_count,
                                RDF, output_file)

            print(f'Job done in {(time() - start)/60:.3f} minutes!')
            print(f'Output file: {output_file}.csv')

        else:
            print('Must choose mono or multi, and select elements to compare.')

    else:
        if args.monocomponent:
            at = args.monocomponent
            print(f'Running 3D RDF between {at} atoms without PBC. Started: {started}.')

            frames_total, Lx, Ly, Lz, nAtTot, nAt, xyz_all = userfile_mono(args.input_file, at)

            nBin, Rcut, RDF = hist_init(Rcut, dr)

            rows = nAtTot + 2

            frames_end = args.frames[1]
            if frames_end == -1:
                frames_end = frames_total

            for frame in range(frames_start, frames_end):
                xyz = xyz_all.iloc[(frame * rows + 2) : ((frame + 1) * rows), :]
                xyz = xyz[xyz['idAt'] == at].to_numpy()
                sample_off_mono(xyz, dr, Rcut, RDF, nAt)
                frames_count += 1

            output_file = args.output_file
            if output_file == None:
                output_file = f'RDF3D_{at}-{at}_dr-{dr:.1f}_Rcut-{Rcut:.1f}'

            normalize_off(dr, Rcut, nBin, frames_count, RDF, output_file)

            print(f'Job done in {(time() - start)/60:.1f} minutes!')
            print(f'Output file: {output_file}.csv')

        elif args.multicomponents:
            at1 = args.multicomponents[0]
            at2 = args.multicomponents[1]
            print(f'Running 3D RDF between {at1} and {at2} atoms without PBC. Started: {started}.')

            frames_total, Lx, Ly, Lz, nAtTot, nAt1, nAt2, xyz_all = userfile_multi(args.input_file, at1, at2)

            nBin, Rcut, RDF = hist_init(Rcut, dr)

            rows = nAtTot + 2

            frames_end = args.frames[1]
            if frames_end == -1:
                frames_end = frames_total

            for frame in range(frames_start, frames_end):
                xyz = xyz_all.iloc[(frame * rows + 2) : ((frame + 1) * rows), :]
                xyz1 = xyz[xyz['idAt'] == at1].to_numpy()
                xyz2 = xyz[xyz['idAt'] == at2].to_numpy()
                sample_off_multi(xyz1, xyz2, dr, Rcut, RDF, nAt1, nAt2)
                frames_count += 1

            output_file = args.output_file
            if output_file == None:
                output_file = f'RDF3D_{at1}-{at2}_dr-{dr:.1f}_Rcut-{Rcut:.1f}'

            normalize_off(dr, Rcut, nBin, frames_count, RDF, output_file)

            print(f'Job done in {(time() - start)/60:.1f} minutes!')
            print(f'Output file: {output_file}.csv')

        else:
            print('Must choose mono or multi, and select elements to compare.')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('input_file', help = "Path to the xsf input file.")

    parser.add_argument('Rcut', type = float, help = "Maximum radius to be \
                        considered.")

    parser.add_argument('dr', type = float, help = "Increment to be considered.")

    parser.add_argument('-f', '--frames', type = int, nargs = 2, default = [0, -1],
                        help = "Choose starting and ending frames to compute.")

    parser.add_argument('-pbc', '--periodic_boundary_conditions',
                        action = 'store_true', help = "Set PBC on.")

    parser.add_argument('-mono', '--monocomponent', help = "Comparison between \
                        the same kind of atom. One argument needed.")

    parser.add_argument('-multi', '--multicomponents', nargs = 2,
                        help = "Comparison between different species. Two \
                        arguments needed.")

    parser.add_argument('-o', '--output_file', help = "Path to the output file. \
                        If not given, the default name will be used.")

    args = parser.parse_args()

    main()
