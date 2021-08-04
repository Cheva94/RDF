#!/usr/bin/env python3.9

'''
    Description: Plots altogether in one plot or one plot each argument.
    Written by: Ignacio J. Chevallier-Boutell.
    Dated: August, 2021.
'''

import argparse
from pandas import read_csv
from cycler import cycler
import matplotlib.pyplot as plt

# plt.rc('font', weight = 'bold', size = 20)
# plt.rc('lines', linewidth = 3)
plt.rc('font', weight = 'bold', size = 15)
V, N, R, Az, = '#00a189', '#fa6200', '#ed3b3b', '#5ca2f7'
plt.rc('axes', labelweight = 'bold', linewidth = 3, prop_cycle = cycler(color=[N, V, 'm', 'y', 'k', 'c', R, Az]))
plt.rc('xtick.major', size = 10, width = 3)
plt.rc('ytick.major', size = 10, width = 3)
plt.rc('legend', loc='best', frameon = True, fancybox = True, fontsize = 15)
plt.rc('figure', autolayout = True)

def main():

    left = args.x_axis[0]
    right = args.x_axis[1]
    bottom = args.y_axis[0]
    top = args.y_axis[1]
    lx = args.x_label
    ly = args.y_label
    hl = args.horizontal_line
    vl = args.vertical_line

    if args.different:

        list = args.different

        for file in list:
            fig, ax = plt.subplots(figsize=(12,10))

            data = read_csv(f'{file}').to_numpy()

            file = file.split('.dat', 1)[0].split('_', 1)[1]

            ax.plot(data[:,0], data[:,1], label = file)
            ax.legend()
            if lx != None:
                ax.set_xlabel(f'{lx}')
            if ly != None:
                ax.set_ylabel(f'{ly}')

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

            plt.savefig(f'{file}.png')
            print(f'Image file: {file}.png')

        print('Job done!')

    elif args.same:
        list = args.same

        fig, ax = plt.subplots(figsize=(12,10))
        for file in list:
            data = read_csv(f'{file}').to_numpy()

            file = file.split('.dat', 1)[0].split('_', 1)[1]

            ax.plot(data[:,0], data[:,1], label = file)
            ax.legend()
            if lx != None:
                ax.set_xlabel(f'{lx}')
            if ly != None:
                ax.set_ylabel(f'{ly}')
            if (left != None) and (right != None):
                ax.set_xlim(float(left), float(right))
            if (bottom != None) and (top != None):
                ax.set_ylim(float(bottom), float(top))

        plt.savefig('same.png')

        print(f'Image file: same.png')

        print('Job done!')

    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('-d', '--different', nargs = '+', help = "Creates one \
                        plot per argument. If just one argument is needed, use \
                        this option.")

    parser.add_argument('-s', '--same', nargs = '+', help = "Creates one plot \
                        overlapping all the arguments.")

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

    args = parser.parse_args()

    main()
