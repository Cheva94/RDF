#!python3.9

import pandas as pd

file_in = 'example31_Si-Pt_VMD_PBC-on.csv'
vmd = pd.read_csv(file_in, header = None, delim_whitespace = True).iloc[:,1].to_numpy()
print(vmd)

with open('on.dat', 'w') as f:
    for row in range(len(vmd)):
        f.write(f'{vmd[row]} \n')
