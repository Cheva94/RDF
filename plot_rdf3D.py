#!/usr/bin/python3.9

'''
    Description: Plots 3D RDF curves.
    Written by: Ignacio J. Chevallier-Boutell.
    Dated: April, 2022.
'''

import argparse
from pandas import read_csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
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
            name = file.split('.csv')[0].split('_', 1)[1]

            peaks, _ = find_peaks(data[:,1], height=0.1, distance=8)
            peaksx, peaksy = data[:,0][peaks], data[:,1][peaks]

            valls, _ = find_peaks(-data[:,1], distance=10)
            vallsx, vallsy = data[:,0][valls], data[:,1][valls]

            fig, ax = plt.subplots()
            ax.plot(data[:,0], data[:,1], label = name)
            ax.plot(peaksx, 1.04*peaksy, lw = 0, marker=11, color='black')
            for i in range(len(peaksx)):
                ax.annotate(f'{peaksx[i]:.2f}', xy = (peaksx[i], 1.05*peaksy[i]), fontsize=10, ha='center')
            ax.plot(vallsx, 0.96*vallsy, lw = 0, marker=10, color='red')
            for i in range(len(vallsx)):
                ax.annotate(f'{vallsx[i]:.2f}', xy = (vallsx[i], 0.88*vallsy[i]), fontsize=10, ha='center', color='red')
            ax.set_xlabel(f'Distance [A]')
            ax.set_ylabel('RDF [3D]')
            ax.legend()

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
        for file in L:
            data = read_csv(f'{file}').to_numpy()
            name = file.split('.csv')[0].split('_', 1)[1]
            F.append(name)

            ax.plot(data[:,0], data[:,1], label = name)
            ax.set_xlabel(f'Distance [A]')
            ax.set_ylabel('RDF [3D]')
            ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1), ncol = 4)

            if (left != None) and (right != None):
                ax.set_xlim(float(left), float(right))
            if (bottom != None) and (top != None):
                ax.set_ylim(float(bottom), float(top))

        ax.xaxis.set_major_locator(MultipleLocator(1))
        ax.xaxis.set_minor_locator(MultipleLocator(0.2))

        plt.savefig(f"RDF3D_{'_'.join(F)}.png")
        print(f"Image file: RDF3D_{'_'.join(F)}.png")

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
