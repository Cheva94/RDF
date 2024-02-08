#!/usr/bin/python3.10

import argparse
from pandas import read_csv
from numpy import sqrt, mean, std, array, argsort

def PBC(dist, length):
    '''
    Correct a distance using Periodic Boundary Conditions (PBC).
    '''

    return dist - length * int(2*dist/length)

def main():

    name = args.input_file
    at1 = args.atoms[0]
    at2 = args.atoms[1]
    Rcut = args.Rcut
    Rcut2 = Rcut * Rcut
    Rmin = args.Rmin
    if Rmin == None:
        Rmin = 0.01
    Rmin2 = Rmin * Rmin

    xsf = read_csv(name, header = None, delim_whitespace = True,
                    names=['idAt', 'rx', 'ry', 'rz', 'fx', 'fy', 'fz'])
    Lx, Ly, Lz = float(xsf.iloc[3,0]), xsf.iloc[4,1], xsf.iloc[5,2]
    xyz = xsf.iloc[8:,0:4].reset_index(drop=True)

    ID = []
    POS = []
    BL = []
    if at1 != at2:
        nAt1 = xyz.iloc[:,0].value_counts()[at1]
        nAt2 = xyz.iloc[:,0].value_counts()[at2]

        xyz1 = xyz[xyz['idAt'] == at1]
        xyz2 = xyz[xyz['idAt'] == at2]

        rx1 = xyz1.iloc[:, 1]
        ry1 = xyz1.iloc[:, 2]
        rz1 = xyz1.iloc[:, 3]

        rx2 = xyz2.iloc[:, 1]
        ry2 = xyz2.iloc[:, 2]
        rz2 = xyz2.iloc[:, 3]

        for i in range(nAt1):
            for j in range(nAt2):
                dx = rx1.iloc[i] - rx2.iloc[j]
                dy = ry1.iloc[i] - ry2.iloc[j]
                dz = rz1.iloc[i] - rz2.iloc[j]

                dx = PBC(dx, Lx)
                dy = PBC(dy, Ly)
                dz = PBC(dz, Lz)

                d2 = dx**2 + dy**2 + dz**2
                if ((d2>= Rmin2) and (d2<= Rcut2)):
                    ID.append(f'-{at1}{rx1.index[i]+1}-{at2}{rx2.index[j]+1}-')
                    POS.append(f'-({rx1.iloc[i]};{ry1.iloc[i]};{rz1.iloc[i]})-({rx2.iloc[j]};{ry2.iloc[j]};{rz2.iloc[j]})-')
                    BL.append(sqrt(d2))

    else:
        nAt = xyz.iloc[:,0].value_counts()[at1]
        xyz = xyz[xyz['idAt'] == at1]
        rx = xyz.iloc[:, 1]
        ry = xyz.iloc[:, 2]
        rz = xyz.iloc[:, 3]

        for i in range(nAt):
            for j in range(i+1, nAt):
                dx = rx.iloc[i] - rx.iloc[j]
                dy = ry.iloc[i] - ry.iloc[j]
                dz = rz.iloc[i] - rz.iloc[j]

                dx = PBC(dx, Lx)
                dy = PBC(dy, Ly)
                dz = PBC(dz, Lz)

                d2 = dx**2 + dy**2 + dz**2
                if ((d2>= Rmin2) and (d2<= Rcut2)):
                    ID.append(f'-{at1}{rx.index[i]+1}-{at2}{rx.index[j]+1}-')
                    POS.append(f'-({rx.iloc[i]};{ry.iloc[i]};{rz.iloc[i]})-({rx.iloc[j]};{ry.iloc[j]};{rz.iloc[j]})-')
                    BL.append(sqrt(d2))

    Summary = f'\nSummary\n  -{at1}-{at2}- = ({mean(BL):.4f} +- {std(BL):.4f}) A\n  Count = {len(BL)}'

    ID = array(ID)
    BL = array(BL)
    POS = array(POS)
    AS = argsort(BL)
    ID = ID[AS]
    BL = BL[AS]
    POS = POS[AS]

    with open(f'BL_{at1}-{at2}.csv', 'w') as f:
        f.write('==== Bond length in Angstroms ==== \n\n')
        f.write('Atoms ID\tDistance\tAtoms Position \n')
        for i in range(len(BL)):
            f.write(f'{ID[i]}\t{BL[i]:.4f}\t{POS[i]} \n')
        f.write(f'{Summary}')

    print(f'{Summary}')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('input_file', help = "Path to the xyz input file.")

    parser.add_argument('Rcut', type = float, help = "Maximum distance to be \
                        considered as bond length.")

    parser.add_argument('atoms', nargs = 2, help = "Atoms to be analyzed.")

    parser.add_argument('--Rmin', type = float,
                        help = "Minimum distance to be considered as bond length.")

    args = parser.parse_args()

    main()
