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
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as mpatches

Verde = '#08a189'
Blanco = '#FFFFFF'
colorGrad = LinearSegmentedColormap.from_list("mycmap", [Blanco, Verde])

def main():

    L = args.input_file
    dxy = args.dxy

    nBinX = int(10.56/dxy) + 1
    xRange = [(binIdxX + 0.5) * dxy for binIdxX in range(nBinX)]
    xLen = len(xRange)

    nBinY = int(14.94/dxy) + 1
    yRange = [(binIdxY + 0.5) * dxy for binIdxY in range(nBinY)]
    yLen = len(yRange)

    for file in L:
        data = read_csv(f'{file}', header=None).to_numpy()
        name = file.split('.csv')[0].split('_')[2]

        _, ax = plt.subplots()

        cresta = ax.axhspan(0, 1.9, color='k', alpha=0.10)
        ax.axhspan(7.5, 9.4, color='k', alpha=0.10)

        mapa = data[:,2].reshape((xLen, yLen)).T
        mapa = np.where(mapa > 0.15*np.max(mapa), mapa, 0.0)
        ax.contour(xRange, yRange, mapa, 40, cmap='summer_r')
        atomo = mpatches.Patch(color=Verde)

        ax.set_xlabel('x [A]')
        ax.set_ylabel('y [A]')
        ax.legend([cresta, atomo], ['Crestas', name], loc='lower center', bbox_to_anchor=(0.5, 1), ncol = 4)

        ax.set_xlim(0, 10.56)
        ax.set_ylim(0, 14.94)

        ax.xaxis.set_major_locator(MultipleLocator(1))
        ax.xaxis.set_minor_locator(MultipleLocator(0.2))

        plt.savefig(f"{file.split('.csv')[0]}.png")
        print(f"Image file: {file.split('.csv')[0]}.png")

    print('Job done!')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('input_file', nargs = '+', help = "Path to the csv files.")

    parser.add_argument('dxy', type = float, help = "Increment to be considered \
                        along x and y axis.")

    args = parser.parse_args()

    main()
