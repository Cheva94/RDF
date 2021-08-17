#!/usr/local/bin/python3.9

'''
    Description: Pre-processing program to consider just the header and the
                positions. By applying PBC, it makes that every particle remains
                inside the given cell.
    Written by: Ignacio J. Chevallier-Boutell.
    Dated: August, 2021.
'''

import argparse
from time import time
from pandas import read_csv
from numpy import array

def main():
    start = time()

    F = args.input_file

    print(f'Running all-in for {F} file.')

    output_file = args.output_file
    if output_file == None:
        output_file = f"{F.split('.xsf', 1)[0]}_all-in"

    with open(f'{F}','r') as orig, open(f'{output_file}.xsf','w') as dest:
        for line in orig.readlines()[:6]:
                dest.write(line)

    xsf = read_csv(F, header = None, delim_whitespace = True,
                    names=['idAt', 'rx', 'ry', 'rz', 'fx', 'fy', 'fz'])

    frames_total = int(xsf.iloc[0,1])
    Lx, Ly, Lz = float(xsf.iloc[3,0]), xsf.iloc[4,1], xsf.iloc[5,2]
    nAtTot = int(xsf.iloc[7,0])
    xyz_all = xsf.iloc[6:,0:4].reset_index(drop=True)
    rows = nAtTot + 2
    L = []

    for frame in range(frames_total):
        xyz = xyz_all.iloc[(frame * rows + 2) : ((frame + 1) * rows), :].to_numpy()
        at = array(xyz[:,0])
        rx = array(xyz[:,1])
        ry = array(xyz[:,2])
        rz = array(xyz[:,3])

        for i in range(nAtTot):

            if rx[i] < 0:
                rx[i] += Lx
            elif rx[i] > Lx:
                rx[i] -= Lx

            if ry[i] < 0:
                ry[i] += Ly
            elif ry[i] > Ly:
                ry[i] -= Ly

            if rz[i] < 0:
                rz[i] += Lz
            elif rz[i] > Lz:
                rz[i] -= Lz

            L.append((at[i], rx[i], ry[i], rz[i]))

    A = array(L)

    with open(f'{output_file}.xsf','a') as f:
        for frame in range(frames_total):
            f.write(f' PRIMCOORD {frame+1} \n')
            f.write(f'{nAtTot} 1 \n')
            xyz = A[(frame * nAtTot) : ((frame + 1) * nAtTot), :]

            for line in range(nAtTot):
                f.write(f'{xyz[line,0]} \t {float(xyz[line,1]):.6f} \t {float(xyz[line,2]):.6f} \t {float(xyz[line,3]):.6f} \n')

    print(f'Job done in {(time() - start):.3f} seconds!')
    print(f'Output file: {output_file}.xsf')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('input_file', help = "Path to the xsf input file.")

    parser.add_argument('-o', '--output_file', help = "Path to the output file. \
                        If not given, the default name will be used.")

    args = parser.parse_args()

    main()
