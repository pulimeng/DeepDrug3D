import argparse
import h5py

import numpy as np

from build_grid import Grid3DBuilder

def site_voxelization(site, r, N, shape):
    """
    Convert the binding site information to numpy array
    """
    site = np.array(site, dtype=np.float64)
    voxel_length = N+1
    voxel_start = -r
    voxel_end = r+1
    coords = site[:,0:3]
    if shape == False:
        print('DFIRE potential included in the voxel representation')
        potentials = site[:,3:]
        voxel = np.zeros(shape=(potentials.shape[1], voxel_length, voxel_length, voxel_length),
            dtype = np.float64)
        cnt = 0
        for x in range(voxel_start, voxel_end+1, 1):
            for y in range(voxel_start, voxel_end+1, 1):
                for z in range(voxel_start, voxel_end+1, 1):
                    temp_voxloc = [x,y,z]
                    distances = np.linalg.norm(coords - temp_voxloc, axis = 1)
                    min_dist = np.min(distances)
                    index = np.where(distances == min_dist)
                    if min_dist < 0.01:
                        voxel[:,x - voxel_start,y - voxel_start,z - voxel_start] = potentials[index,:]
                        cnt += 1
                    else:
                        voxel[:,x - voxel_start,y - voxel_start,z - voxel_start] = np.ones((14,))
    else:
        print('Binary occupation only for voxel representation')
        potentials = np.ones((site.shape[0],1))
        voxel = np.zeros(shape=(1, voxel_length, voxel_length, voxel_length),
            dtype = np.float64) 
        cnt = 0
        for x in range(voxel_start, voxel_end+1, 1):
            for y in range(voxel_start, voxel_end+1, 1):
                for z in range(voxel_start, voxel_end+1, 1):
                    temp_voxloc = [x,y,z]
                    distances = np.linalg.norm(coords - temp_voxloc, axis = 1)
                    min_dist = np.min(distances)
                    index = np.where(distances == min_dist)
                    if min_dist < 0.01:
                        voxel[:,x - voxel_start,y - voxel_start,z - voxel_start] = potentials[index,:]
                        cnt += 1
                    else:
                        voxel[:,x - voxel_start,y - voxel_start,z - voxel_start] = np.zeros((1,))
    return voxel

class Vox3DBuilder(object):
    """
    This class convert the pdb file to the voxel representation for the input
    of deep learning architecture. The conversion is around 30 mins.
    """
    @staticmethod
    def voxelization(pdb_path, aux_input_path, r, N, output_folder, shape):
        print('Generating pocket grid representation')
        pocket_grid = Grid3DBuilder.build(pdb_path, aux_input_path, r, N, output_folder, shape)
        print('Converting to numpy array')
        pocket_voxel = site_voxelization(pocket_grid, r, N, shape)
        pocket_voxel = np.expand_dims(pocket_voxel, axis=0)
        with h5py.File('{}/{}.h5'.format(output_folder, pdb_path[:-4]), 'w') as f:
            f.create_dataset('X', data=pocket_voxel, compression='gzip', compression_opts=9)
        return pocket_voxel

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--f', type=str, required=True, help='Input pdb file name')
    parser.add_argument('--a', type=str, required=True, help='Input auxilary file name')
    parser.add_argument('--o', type=str, required=True, help='Output folder name')
    parser.add_argument('--r', type=int, required=True, help='Grid radius')
    parser.add_argument('--n', type=int, required=True, help='Number of points along diameter')
    parser.add_argument('--s', dest='shape', action='store_true', help='Return shape of grid only')
    parser.add_argument('--p', dest='shape', action='store_false')
    parser.set_defaults(shape=False)

    opt = parser.parse_args()
    Vox3DBuilder().voxelization(opt.f, opt.a, opt.r, opt.n, opt.o, opt.shape)