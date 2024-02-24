#!/usr/bin/python3.10

import argparse
from time import time
from pandas import read_csv
from numpy import array

def main():
    start = time()

    F = args.input_file
    Zcoord = args.Zcoord

    print(f'Running centrarPoro for {F} file.')

    output_file = args.output_file
    if output_file == None:
        output_file = f"{F.split('.xsf', 1)[0]}_PoroCentrado"

    with open(f'{F}','r') as orig, open(f'{output_file}.xsf','w') as dest:
        for line in orig.readlines()[:6]:
                dest.write(line)

    xsf = read_csv(F, header = None, delim_whitespace = True,
                    names=['idAt', 'rx', 'ry', 'rz', 'fx', 'fy', 'fz'])

    frames_total = int(xsf.iloc[0,1])
    nAtTot = int(xsf.iloc[7,0])
    rows = nAtTot + 2
    xyz_all = xsf.iloc[6:,0:4].reset_index(drop=True)
    L = []
    oldFrame = 0

    for frame in range(frames_total):
        xyz = xyz_all.iloc[(frame * rows + 2) : ((frame + 1) * rows), :].to_numpy()
        at = array(xyz[:,0])
        rx = array(xyz[:,1]).astype(float)
        ry = array(xyz[:,2]).astype(float)
        rz = array(xyz[:,3]).astype(float)

        for i in range(nAtTot):
            rz[i] += Zcoord
            L.append((at[i], rx[i], ry[i], rz[i]))

        cutmod = 5000
        if (frame+1) % cutmod  == 0:
            print(f'dump: {frame+1}')

            A = array(L)

            with open(f'{output_file}.xsf','a') as f:
                for partframe in range(cutmod):
                    f.write(f' PRIMCOORD {oldFrame+partframe+1} \n')
                    f.write(f'{nAtTot} 1 \n')
                    xyzw = A[(partframe * nAtTot) : ((partframe + 1) * nAtTot), :]

                    for line in range(nAtTot):
                        f.write(f'{xyzw[line,0]} \t {float(xyzw[line,1]):.6f} \t {float(xyzw[line,2]):.6f} \t {float(xyzw[line,3]):.6f} \n')

            L = []
            oldFrame += cutmod

    A = array(L)

    with open(f'{output_file}.xsf','a') as f:
        for lastframe in range(frames_total-oldFrame):
            f.write(f' PRIMCOORD {oldFrame+lastframe+1} \n')
            f.write(f'{nAtTot} 1 \n')
            xyz = A[(lastframe * nAtTot) : ((lastframe + 1) * nAtTot), :]

            for line in range(nAtTot):
                f.write(f'{xyz[line,0]} \t {float(xyz[line,1]):.6f} \t {float(xyz[line,2]):.6f} \t {float(xyz[line,3]):.6f} \n')

    print(f'Job done in {(time() - start)/60:.3f} minutes!')
    print(f'Output file: {output_file}.xsf')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('input_file', help = "Path to the xsf input file.")
    parser.add_argument('Zcoord', help = "Change in Z coordnate.", type=float)

    parser.add_argument('-o', '--output_file', help = "Path to the output file. \
                        If not given, the default name will be used.")

    args = parser.parse_args()

    main()