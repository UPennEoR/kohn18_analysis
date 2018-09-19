import argparse
import pyuvdata
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument(
    '-f',
    '--files',
    help='files to be converted from miriad to hdf5 file format',
    nargs='*',
    required=True)
args = parser.parse_args()
files = np.array(sorted(args.files))

for miriad_file in files:
    uvh5_file = miriad_file + '.uvh5'
    print miriad_file + ' --> ' + uvh5_file
    uvd = pyuvdata.UVData()
    uvd.read_miriad(miriad_file)
    uvd.write_uvh5(uvh5_file)