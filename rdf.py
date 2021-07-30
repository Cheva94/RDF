#!/usr/bin/env python3.9

'''
    Calculation: Radial Distribution Function (RDF).
    Description: using Angstrom.
    Written by: Ignacio J. Chevallier-Boutell.
    Dated: July, 2021.
'''

import argparse
from core_rdf import *

def main():
    if args.periodic_boundary_conditions:
        if args.monocomponent:

            dr = args.dr
            Rcut = args.Rcut
            at = args.monocomponent

            total_frames, Lx, Ly, Lz, nAtTot, nAt, xyz_all = user_file_mono(args.input_file, at)

            Lmin = 0.5 * min(Lx, Ly, Lz)
            if Rcut > Lmin:
                print(f'Cannot choose Rcut greater than {Lmin:.2f}. This will be the new Rcut value.')
                Rcut = Lmin

            nBin, Rcut, RDF = hist_init(dr, Rcut)

            frames_count = 0
            rows = nAtTot + 2
            for frame in range(total_frames):
                xyz = xyz_all.iloc[(frame*rows + 2):((frame+1)*rows) , :]
                xyz = xyz[xyz['idAt'] == at].to_numpy()
                mono_on_sample(Lx, Ly, Lz, xyz, dr, Rcut, RDF, nAt)
                frames_count += 1

            output_file = args.output_file
            if output_file == None:
                output_file = f'rdf_{at}-{at}_PBC-on'

            mono_on_normalize(Lx, Ly, Lz, nAt, dr, nBin, frames_count, RDF, output_file)

            print(f'Job done! RDF between {at} and {at} was calculated with PBC. \nOutput file: {output_file}.dat')

        elif args.multicomponents:

            dr = args.dr
            Rcut = args.Rcut
            at1 = args.multicomponents[0]
            at2 = args.multicomponents[1]

            total_frames, Lx, Ly, Lz, nAtTot, nAt1, nAt2, xyz_all = user_file_multi(args.input_file, at1, at2)

            Lmin = 0.5 * min(Lx, Ly, Lz)
            if Rcut > Lmin:
                print(f'Cannot choose Rcut greater than {Lmin:.2f}. This will be the new Rcut value.')
                Rcut = Lmin

            nBin, Rcut, RDF = hist_init(dr, Rcut)

            frames_count = 0
            rows = nAtTot + 2
            for frame in range(total_frames):
                xyz = xyz_all.iloc[(frame*rows + 2):((frame+1)*rows) , :]
                xyz1 = xyz[xyz['idAt'] == at1].to_numpy()
                xyz2 = xyz[xyz['idAt'] == at2].to_numpy()
                multi_on_sample(Lx, Ly, Lz, xyz1, xyz2, dr, Rcut, RDF, nAt1, nAt2)
                frames_count += 1

            output_file = args.output_file
            if output_file == None:
                output_file = f'rdf_{at1}-{at2}_PBC-on'

            multi_on_normalize(Lx, Ly, Lz, nAt1, nAt2, dr, nBin, frames_count, RDF, output_file)

            print(f'Job done! RDF between {at1} and {at2} was calculated with PBC. \nOutput file: {output_file}.dat')

        else:
            print('Must choose mono or multi, and select elements to compare.')

    else:
        if args.monocomponent:

            dr = args.dr
            Rcut = args.Rcut
            at = args.monocomponent

            total_frames, Lx, Ly, Lz, nAtTot, nAt, xyz_all = user_file_mono(args.input_file, at)

            nBin, Rcut, RDF = hist_init(dr, Rcut)

            frames_count = 0
            rows = nAtTot + 2
            for frame in range(total_frames):
                xyz = xyz_all.iloc[(frame*rows + 2):((frame+1)*rows) , :]
                xyz = xyz[xyz['idAt'] == at].to_numpy()
                mono_off_sample(xyz, dr, Rcut, RDF, nAt)
                frames_count += 1

            output_file = args.output_file
            if output_file == None:
                output_file = f'rdf_{at}-{at}_PBC-off'

            mono_off_normalize(nAt, dr, nBin, frames_count, RDF, output_file)

            print(f'Job done! RDF between {at} and {at} was calculated without PBC. \nOutput file: {output_file}.dat')

        elif args.multicomponents:

            dr = args.dr
            Rcut = args.Rcut
            at1 = args.multicomponents[0]
            at2 = args.multicomponents[1]

            total_frames, Lx, Ly, Lz, nAtTot, nAt1, nAt2, xyz_all = user_file_multi(args.input_file, at1, at2)

            nBin, Rcut, RDF = hist_init(dr, Rcut)

            frames_count = 0
            rows = nAtTot + 2
            for frame in range(total_frames):
                xyz = xyz_all.iloc[(frame*rows + 2):((frame+1)*rows) , :]
                xyz1 = xyz[xyz['idAt'] == at1].to_numpy()
                xyz2 = xyz[xyz['idAt'] == at2].to_numpy()
                multi_off_sample(xyz1, xyz2, dr, Rcut, RDF, nAt1, nAt2)
                frames_count += 1

            output_file = args.output_file
            if output_file == None:
                output_file = f'rdf_{at1}-{at2}_PBC-off'

            multi_off_normalize(nAt1, nAt2, dr, nBin, frames_count, RDF, output_file)

            print(f'Job done! RDF between {at1} and {at2} was calculated without PBC. \nOutput file: {output_file}.dat')

        else:
            print('Must choose mono or multi, and select elements to compare.')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('-pbc', '--periodic_boundary_conditions', action = 'store_true', help = "quiere PBC")

    parser.add_argument('-mono', '--monocomponent', help = "tira mono o multi")

    parser.add_argument('-multi', '--multicomponents', nargs = 2, help = "tira mono o multi")

    parser.add_argument('dr', type = float, help = "increment")

    parser.add_argument('Rcut', type = float, help = "maximum radius")

    parser.add_argument('input_file', help = "path del xsf")

    parser.add_argument('-o', '--output_file', help = "path del output")

    args = parser.parse_args()

    main()
