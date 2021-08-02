#!/usr/bin/env python3.9

'''
    Description: Plots altogether in one plot or one plot each argument.
    Written by: Ignacio J. Chevallier-Boutell.
    Dated: July, 2021.
'''

import argparse
from pandas import read_csv
from cycler import cycler
import matplotlib.pyplot as plt

plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams.update({'font.size': 35})

plt.rcParams['lines.linewidth'] = 3
plt.rcParams['axes.linewidth'] = 3
plt.rcParams['xtick.major.size'] = 10
plt.rcParams['xtick.major.width'] = 3
plt.rcParams['ytick.major.size'] = 10
plt.rcParams['ytick.major.width'] = 3

V = '#00a189'
N = '#fa6200'
R = '#ed3b3b'
Az = '#5ca2f7'
Az_claro = '#97C4FA'
Am = '#ebe842'

plt.rcParams['axes.prop_cycle'] = cycler(color=[V, N, R, Az, Am, Az_claro])

def main():
    if args.different:
        left = args.x_axis[0]
        right = args.x_axis[1]
        bottom = args.y_axis[0]
        top = args.y_axis[1]

        list = args.different

        for file in list:
            fig, ax = plt.subplots(figsize=(12,10))

            data = read_csv(f'{file}').to_numpy()

            file = file.split('.', 1)[0]
            ylabel = file.split('_',1)[0]

            ax.plot(data[:,0], data[:,1], label = file)
            ax.set_ylabel(f'{ylabel}')
            ax.set_xlabel('Distance [A]')
            ax.legend(loc='upper right', fontsize=15)

            if (left != None) and (right != None):
                ax.set_xlim(float(left), float(right))

            if (bottom != None) and (top != None):
                ax.set_ylim(float(bottom), float(top))

            plt.tight_layout()

            plt.savefig(f'{file}.png')
            print(f'Image file: {file}.png')

        print('Job done!')

    elif args.same:
        left = args.x_axis[0]
        right = args.x_axis[1]
        bottom = args.y_axis[0]
        top = args.y_axis[1]

        list = args.same

        fig, ax = plt.subplots(figsize=(12,10))
        for file in list:
            data = read_csv(f'{file}').to_numpy()

            file = file.split('.', 1)[0]
            ylabel = file.split('_',1)[0]

            ax.plot(data[:,0], data[:,1], label = file)
            ax.set_ylabel(f'{ylabel}')
            ax.set_xlabel('Distance [A]')
            ax.legend(loc='upper right', fontsize=15)

            if (left != None) and (right != None):
                ax.set_xlim(float(left), float(right))

            if (bottom != None) and (top != None):
                ax.set_ylim(float(bottom), float(top))

            plt.tight_layout()

            plt.savefig(f'{ylabel}.png')

        print(f'Image file: {ylabel}.png')

        print('Job done!')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('-d', '--different', nargs = '+', help = "Creates one \
                        plot per argument. If just one argument is needed, use \
                        this option.")

    parser.add_argument('-s', '--same', nargs = '+', help = "Creates one plot \
                        overlapping all the arguments.")

    parser.add_argument('-x', '--x_axis', nargs = 2, default = [None, None],
                        help = "Choose range for X axis .")

    parser.add_argument('-y', '--y_axis', nargs = 2, default = [None, None],
                        help = "Choose range for Y axis .")

    args = parser.parse_args()

    main()
