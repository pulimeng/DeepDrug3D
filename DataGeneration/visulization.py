import argparse

import numpy as np
import pandas as pd
import string

def makeLine(coord, potE, k, num_chr_dict):
    atom = 'ATOM'
    atom_sn = str(k+1)
    atom_name = 'D1'
    hundred, one = np.divmod(k,676)
    ten, one = np.divmod(k-hundred*676,26)
    res_name = num_chr_dict[hundred] + num_chr_dict[ten] + num_chr_dict[one]
    E = potE
    x = '{:.3f}'.format(round(coord[0], 3))
    y = '{:.3f}'.format(round(coord[1], 3))
    z = '{:.3f}'.format(round(coord[2], 3))
    EE = '{:.2f}'.format(round(E,2))
    string = atom + ' '*2 + '{:>5}'.format(atom_sn) + ' ' + '{:4}'.format(atom_name) + ' ' \
            + '{:>3}'.format(res_name) + ' '*2 + '   1' + ' '*4 + '{:>8}'.format(x) + '{:>8}'.format(y) + '{:>8}'.format(z) \
            + '{:>6}'.format('1.00') + '{:>6}'.format(EE) + ' '*8 + '\n'
    return string

def main(opt):
    num_range = np.linspace(0,25,26, dtype = int)
    chr_range = list(string.ascii_uppercase[:27])
    num_chr_dict = dict(zip(num_range,chr_range))
    
    site = pd.read_csv(opt.i)
    site = site.to_numpy()
    with open ('{}_grid_channel_{}.pdb'.format(opt.i[:-5], opt.c),'w') as in_strm:
        for k in range(len(site)):
            temp_coords = site[k,:]
            cl = temp_coords[opt.c]
            temp_string = makeLine(temp_coords, cl, k, num_chr_dict)
            in_strm.write(temp_string)
    print('Number of points is ' + str(len(site)))
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--i', type=str, help='Input grid file')
    parser.add_argument('--c', type=int, help='Channel to visualize')
    opt = parser.parse_args()
    main(opt)