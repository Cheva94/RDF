#!/usr/local/bin/python3.9

'''
    Description: Plots altogether in one plot.
    Written by: Ignacio J. Chevallier-Boutell.
    Dated: August, 2021.
'''

import argparse
from pandas import read_csv
import matplotlib.pyplot as plt

def main():

    plt.style.use('strucan')

    left = args.x_axis[0]
    right = args.x_axis[1]
    bottom = args.y_axis[0]
    top = args.y_axis[1]
    lx = args.x_label
    ly = args.y_label
    hl = args.horizontal_line
    vl = args.vertical_line

    L = args.input_file
    F = []

    fig, ax = plt.subplots()
    for file in L:

        data = read_csv(f'{file}').to_numpy()

        name = file.split('.csv')[0].split('_', 1)[1]

        if args.pdf == False:
            if args.rdf3d:
                name = name.split('_')
                F.append(name[0])
                name = f"{name[0]} dr:{name[2].split('-')[1]}"
            elif args.rdf2d:
                name = name.split('_')
                F.append(name[0])
                name = f"{name[0]} dr:{name[3].split('-')[1]} H:{name[1].split('-', 1)[1]}"
            elif args.hdf:
                F.append(name)
            else:
                print('Must choose between --rdf3d (-R3), --rdf2d (-R2) and --hdf (-H).')
                exit()

            ax.plot(data[:,0], data[:,1], label = name)
        else:
            F.append(name)
            
            ax.scatter(data[:,0], data[:,1], label = name, s = 10)

    cols = len(L)

    if (cols % 2) == 0:
       cols = int(cols*0.5)
    else:
       cols = int(cols*0.5) + 1

    ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1), ncol = cols)

    if (left != None) and (right != None):
        ax.set_xlim(float(left), float(right))
    if (bottom != None) and (top != None):
        ax.set_ylim(float(bottom), float(top))

    if vl != None:
        for line in vl:
            ax.axvline(float(line), c = 'black', lw = 1.5, ls = ':')
    if hl != None:
        for line in hl:
            ax.axhline(float(line), c = 'black', lw = 1.5, ls = ':')

    if args.rdf3d:

        if lx == None:
            ax.set_xlabel(f'Distance [A]')
        else:
            ax.set_xlabel(f'{lx}')

        if ly == None:
            ax.set_ylabel(f'g(r) [3D]')
        else:
            ax.set_ylabel(f'{ly}')

        plt.savefig(f"RDF3D_{'_'.join(F)}.png")
        print(f"Image file: RDF3D_{'_'.join(F)}.png")

    elif args.rdf2d:

        if lx == None:
            ax.set_xlabel(f'Distance [A]')
        else:
            ax.set_xlabel(f'{lx}')

        if ly == None:
            ax.set_ylabel(f'g(r) [2D]')
        else:
            ax.set_ylabel(f'{ly}')

        plt.savefig(f"RDF2D_{'_'.join(F)}.png")
        print(f"Image file: RDF2D_{'_'.join(F)}.png")

    elif args.hdf:

        if lx == None:
            ax.set_xlabel(f'Height [A]')
        else:
            ax.set_xlabel(f'{lx}')

        if ly == None:
            ax.set_ylabel(f'h(z)')
        else:
            ax.set_ylabel(f'{ly}')

        plt.savefig(f"plt-s_HDF-{'-'.join(F)}.png")
        print(f"Image file: plt-s_HDF-{'-'.join(F)}.png")

    else:

        if lx == None:
            ax.set_xlabel(f'x [A]')
        else:
            ax.set_xlabel(f'{lx}')

        if ly == None:
            ax.set_ylabel(f'y [A]')
        else:
            ax.set_ylabel(f'{ly}')

        plt.savefig(f"PDF_{'_'.join(F)}.png")
        print(f"Image file: PDF_{'_'.join(F)}.png")

    print('Job done!')

    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('input_file', nargs = '+', help = "Path to the csv files.")

    parser.add_argument('-x', '--x_axis', nargs = 2, default = [None, None],
                        help = "Choose range for X axis.")

    parser.add_argument('-y', '--y_axis', nargs = 2, default = [None, None],
                        help = "Choose range for Y axis.")

    parser.add_argument('-lx', '--x_label', default = None,
                        help = "Choose label for X axis.")

    parser.add_argument('-ly', '--y_label', default = None,
                        help = "Choose label for Y axis.")

    parser.add_argument('-hl', '--horizontal_line', nargs = '+', default = None,
                        help = "Add horizontal line.")

    parser.add_argument('-vl', '--vertical_line', nargs = '+', default = None,
                        help = "Add vertical line.")

    parser.add_argument('-R3', '--rdf3d', action = 'store_true', help = "Creates one plot \
                        per argument for 3D RDF files.")

    parser.add_argument('-R2', '--rdf2d', action = 'store_true', help = "Creates one plot \
                        per argument for 2D RDF files.")

    parser.add_argument('-H', '--hdf', action = 'store_true', help = "Creates one plot \
                        per argument for HDF files.")

    parser.add_argument('-P', '--pdf', action = 'store_true', help = "Creates one plot \
                        per argument for PDF files.")

    args = parser.parse_args()

    main()
