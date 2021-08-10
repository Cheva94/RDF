 #!/usr/bin/env python3.9

'''
    Calculation: Plane Distribution Function (PDF).
    Description: Determines de PDF over XY plane when given a xsf file.
                Everything is in Angstrom.
    Written by: Ignacio J. Chevallier-Boutell.
    Dated: August, 2021.
'''

import argparse
from core import user_file_mono, hist_init2d
from numpy import array
from time import time

def main():
    start = time() # starting wall time

    at = args.at
    dxy = args.dxy
    # dh = args.dh
    Hmin = args.Hcut[0]
    Hmax = args.Hcut[1]
    print(f'Running PDF for {at}.')

    total_frames, Lx, Ly, Lz, nAtTot, nAt, xyz_all = user_file_mono(args.input_file, at)

    nBinX, nBinY, Lx, Ly, PDF = hist_init2d(0, Lx, 0, Ly, dxy)

    frames_count = 0
    rows = nAtTot + 2

    frame_start = args.frames[0]
    if frame_start != 0:
        frame_start -= 1
    frame_end = args.frames[1]
    if frame_end == -1:
        frame_end = total_frames

    for frame in range(frame_start, frame_end):
        xyz = xyz_all.iloc[(frame*rows + 2):((frame+1)*rows) , :]
        xyz = xyz[(xyz['idAt'] == at) & (Hmin <= xyz['rz']) & (xyz['rz'] < Hmax)].to_numpy()

        dx = array(xyz[:,1])
        dy = array(xyz[:,2])

        nAt = len(xyz)

        for i in range(nAt):
            if dx[i] <= 0:
                dx[i] += Lx
            elif dx[i] > Lx:
                dx[i] -= Lx
            binIdxX = int(dx[i]/dxy)

            if dy[i] <= 0:
                dy[i] += Ly
            elif dy[i] > Ly:
                dy[i] -= Ly
            binIdxY = int(dy[i]/dxy)

            PDF[binIdxX, binIdxY] += 1

        frames_count += 1

    output_file = args.output_file
    if output_file == None:
        output_file = f'PDF_{at}'#_dh-{dh}'

    with open(f'{output_file}.csv', 'w') as f:
        for binIdxX in range(nBinX):
            for binIdxY in range(nBinY):
                x = (binIdxX + 0.5) * dxy
                y = (binIdxY + 0.5) * dxy
                PDF[binIdxX, binIdxY] /= frames_count
                if PDF[binIdxX, binIdxY] == 0:
                    continue
                else:
                    f.write(f'{x:.2f}, {y:.2f}, {PDF[binIdxX, binIdxY]:.4f} \n')

    elapsed = time() - start # elapsed wall time
    print(f'Job done in {elapsed:.3f} seconds!')
    print(f'Output file: {output_file}.csv')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('input_file', help = "Path to the xsf input file.")

    parser.add_argument('at', help = "Atom to be analyzed.")

    parser.add_argument('dxy', type = float, help = "Increment to be considered \
                        along x and y axis.")

    parser.add_argument('Hcut', type = float, nargs = 2, default = [0, -1],
                        help = "Minimum and maximum heights to be considered.")

    # parser.add_argument('dh', type = float, help = "Width of the slab to be \
    #                     considered")

    parser.add_argument('-o', '--output_file', help = "Path to the output file. \
                        If not given, the default name will be used.")

    parser.add_argument('-f', '--frames', nargs = 2, default = [0, -1], type = int,
                        help = "Choose starting and ending frames to compute.")

    args = parser.parse_args()

    main()
