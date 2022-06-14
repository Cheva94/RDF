#!/usr/bin/python3.10

'''
    Description: Plots PDF curves.
    Written by: Ignacio J. Chevallier-Boutell.
    Dated: April, 2022.
'''

import argparse
from pandas import read_csv
import numpy as np
import matplotlib.pyplot as plt
# from scipy.signal import find_peaks
from matplotlib.ticker import MultipleLocator

def main():

    left = args.x_axis[0]
    right = args.x_axis[1]
    bottom = args.y_axis[0]
    top = args.y_axis[1]
    L = args.input_file

    if args.multiPlot:

        for file in L:
            data = read_csv(f'{file}').to_numpy()
            name = file.split('.csv')[0].split('_')[1]

            fig, ax = plt.subplots()

            ax.axhspan(0.125, 1.902, color='k', alpha=0.10, label = "Ridges")
            ax.axhspan(7.449, 9.286, color='k', alpha=0.10)

            disp = ax.scatter(data[:,0], data[:,1], c = data[:,2], cmap = 'summer_r', s=10, label = name)
            fig.colorbar(disp)
            ax.set_xlabel('x [A]')
            ax.set_ylabel('y [A]')
            ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1), ncol = 4)

            if (left != None) and (right != None):
                ax.set_xlim(float(left), float(right))
            if (bottom != None) and (top != None):
                ax.set_ylim(float(bottom), float(top))

            ax.xaxis.set_major_locator(MultipleLocator(1))
            ax.xaxis.set_minor_locator(MultipleLocator(0.2))

            plt.savefig(f"{file.split('.csv')[0]}.png")
            print(f"Image file: {file.split('.csv')[0]}.png")

    elif args.onePlot:

        F = []
        fig, ax = plt.subplots()
        
        ax.axhspan(0.125, 1.902, color='k', alpha=0.10, label = "Ridges")
        ax.axhspan(7.449, 9.286, color='k', alpha=0.10)

        for file in L:
            data = read_csv(f'{file}').to_numpy()
            name = file.split('.csv')[0].split('_')[1]
            F.append(name)

            ax.scatter(data[:,0], data[:,1], label = name, s = 10)
            ax.set_xlabel('x [A]')
            ax.set_ylabel('y [A]')
            ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1), ncol = 4)

            if (left != None) and (right != None):
                ax.set_xlim(float(left), float(right))
            if (bottom != None) and (top != None):
                ax.set_ylim(float(bottom), float(top))

        ax.xaxis.set_major_locator(MultipleLocator(1))
        ax.xaxis.set_minor_locator(MultipleLocator(0.2))

        plt.savefig(f"PDF_{'_'.join(F)}.png")
        print(f"Image file: PDF_{'_'.join(F)}.png")

    else:
        print('Must choose between one plot (-one) or multiple plots (-multi).')
        exit()

    print('Job done!')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('input_file', nargs = '+', help = "Path to the csv files.")

    parser.add_argument('-x', '--x_axis', nargs = 2, default = [None, None],
                        help = "Choose range for X axis.")

    parser.add_argument('-y', '--y_axis', nargs = 2, default = [None, None],
                        help = "Choose range for Y axis.")

    parser.add_argument('-one', '--onePlot', action = 'store_true', help = "Plots altogether in one plot.")

    parser.add_argument('-multi', '--multiPlot', action = 'store_true', help = "Plots each argument on one plot.")

    args = parser.parse_args()

    main()
