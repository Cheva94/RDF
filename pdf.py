 #!/usr/bin/env python3.9

'''
<<<<<<< HEAD
    Calculation: Plane Distribution Function (PDF).
    Description: Determines de PDF over XY plane when given a xsf file.
=======
    Calculation: Height Distribution Function (HDF).
    Description: Determines de HDF along Z axis when given a xsf file.
>>>>>>> e6e93410b4e69b1ce0a2119ebe604732160f7da0
                Everything is in Angstrom.
    Written by: Ignacio J. Chevallier-Boutell.
    Dated: August, 2021.
'''

import argparse
<<<<<<< HEAD
from core import user_file_mono, hist_init2d
=======
from core import user_file_mono, hist_init
>>>>>>> e6e93410b4e69b1ce0a2119ebe604732160f7da0
from numpy import array
from time import time

def main():
    start = time() # starting wall time

<<<<<<< HEAD
    at = args.at
    dxy = args.dxy
    # dh = args.dh
    Hmin = args.Hcut[0]
    Hmax = args.Hcut[1]
    print(f'Running PDF for {at}.')

    total_frames, Lx, Ly, Lz, nAtTot, nAt, xyz_all = user_file_mono(args.input_file, at)

    nBinX, nBinY, Lx, Ly, PDF = hist_init2d(0, Lx, 0, Ly, dxy)
=======
    dh = args.dh
    at = args.at
    Hmin = args.Hcut[0]
    Hmax = args.Hcut[1]
    print(f'Running HDF for {at}.')

    total_frames, Lx, Ly, Lz, nAtTot, nAt, xyz_all = user_file_mono(args.input_file, at)

    if Hmax == -1:
        Hmax = Lz

    nBin, Hmax, HDF = hist_init(Hmin, Hmax, dh)
>>>>>>> e6e93410b4e69b1ce0a2119ebe604732160f7da0

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
<<<<<<< HEAD
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
=======
        xyz = xyz[(xyz['idAt'] == at) & (Hmin <= xyz['rz']) & (xyz['rz'] <= Hmax)].to_numpy()

        dz = array(xyz[:,3])

        nAt = len(xyz)
        for i in range(nAt):
            binIdx = int((dz[i] - Hmin)/dh)
            HDF[binIdx] += 1
>>>>>>> e6e93410b4e69b1ce0a2119ebe604732160f7da0

        frames_count += 1

    output_file = args.output_file
    if output_file == None:
<<<<<<< HEAD
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
=======
        output_file = f'HDF_{at}_dh-{dh}'

    with open(f'{output_file}.dat', 'w') as f:
        for binIdx in range(nBin):
            h = (binIdx + 0.5) * dh + Hmin
            HDF[binIdx] /= frames_count
            f.write(f'{h:.2f}, {HDF[binIdx]:.4f} \n')

    elapsed = time() - start # elapsed wall time
    print(f'Job done in {elapsed:.3f} seconds!')
    print(f'Output file: {output_file}.dat')
>>>>>>> e6e93410b4e69b1ce0a2119ebe604732160f7da0

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('input_file', help = "Path to the xsf input file.")

    parser.add_argument('at', help = "Atom to be analyzed.")

<<<<<<< HEAD
    parser.add_argument('dxy', type = float, help = "Increment to be considered \
                        along x and y axis.")

    parser.add_argument('Hcut', type = float, nargs = 2, default = [0, -1],
                        help = "Minimum and maximum heights to be considered.")

    # parser.add_argument('dh', type = float, help = "Width of the slab to be \
    #                     considered")
=======
    parser.add_argument('dh', type = float, help = "Increment to be considered.")
>>>>>>> e6e93410b4e69b1ce0a2119ebe604732160f7da0

    parser.add_argument('-o', '--output_file', help = "Path to the output file. \
                        If not given, the default name will be used.")

<<<<<<< HEAD
    parser.add_argument('-f', '--frames', nargs = 2, default = [0, -1], type = int,
                        help = "Choose starting and ending frames to compute.")

=======
    parser.add_argument('-f', '--frames', type = int, nargs = 2, default = [0, -1],
                        help = "Choose starting and ending frames to compute.")

    parser.add_argument('-H', '--Hcut', type = float, nargs = 2, default = [0, -1],
                        help = "Minimum and maximum heights to be considered.")

>>>>>>> e6e93410b4e69b1ce0a2119ebe604732160f7da0
    args = parser.parse_args()

    main()
