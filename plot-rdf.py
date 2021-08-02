#!/usr/bin/env python3.9

'''
    Description: Plots RDF. It can be altogether in one plot or one plot each RDF.
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
        list = args.different

        for file in list:
            fig, ax = plt.subplots(figsize=(12,10))

            data = read_csv(f'{file}').to_numpy()

            file = file.split('.', 1)[0]
            title = f"{file.split('_')[1]}    {file.split('_')[2].split('-')[0]}: {file.split('_')[2].split('-')[1]}"

            ax.set_title(f'{title}')
            ax.plot(data[:,0], data[:,1])
            ax.set_ylabel('g(r)')
            ax.set_xlabel('r [A]')

            plt.tight_layout()

            plt.savefig(f'{file}.png')
            print(f'Image file: {file}.png')

        print('Job done!')

    elif args.same:
        list = args.same

        fig, ax = plt.subplots(figsize=(12,10))
        for file in list:
            data = read_csv(f'{file}').to_numpy()

            file = file.split('.', 1)[0].split('_',1)[1].split('_')
            file = f"{file[0]}    {file[1].split('-')[0]}: {file[1].split('-')[1]}"

            ax.plot(data[:,0], data[:,1], label = file)
            ax.set_ylabel('g(r)')
            ax.set_xlabel('r [A]')
            ax.legend(loc='upper right', fontsize=15)

            plt.tight_layout()

            plt.savefig(f'rdf.png')

        print(f'Image file: rdf.png')

        print('Job done!')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('-d', '--different', nargs = '+', help = "Creates one \
                        plot per argument. If just one argument is needed, use \
                        this option.")

    parser.add_argument('-s', '--same', nargs = '+', help = "Creates one plot \
                        overlapping all the arguments.")

    args = parser.parse_args()

    main()
