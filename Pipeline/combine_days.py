"""This script will combine all given uvf5 files into one uvf5 file; expects they are all of one pol.
This was used to create one file for each pol that has info for all 8 days in it, used after the
appropriate LSTs had already been selected from each day and gathered into 1 hdf5 file for that day."""

import os
import argparse
import numpy as np
import astropy.units as u
import astropy.coordinates as aco
from pyuvdata import UVData

parser = argparse.ArgumentParser()
parser.add_argument(
    '-F',
    '--files',
    help='Designate the hdf5 files to be concatenated in time.',
    nargs='*',
    required=True)
parser.add_argument(
    '-S',
    '--savepath',
    help='Designate the path where the new hdf5 files will be saved. Default is path to data files.')
args = parser.parse_args()

#Formatting command line arguments
files = np.array(sorted(args.files))
if args.savepath is None:
    savepath = os.path.dirname(args.files[0])
else:
    savepath = args.savepath
    os.system('mkdir -p {}'.format(savepath))

#Load in and combine the files
uvd = UVData()
print 'Reading in: ' + files[0]
uvd.read_uvh5(
    files[0],
    ant_str='cross')
for dfile in files[1:]:
    print 'Reading in: ' + dfile
    uvdi = UVData()
    uvdi.read_uvh5(
        dfile, 
        ant_str='cross')
    uvd += uvdi
    
#Name and save the new object
JDi = files[0].split('/')[-1].split('.')[1]
JDf = files[-1].split('/')[-1].split('.')[1]
name_parts = files[0].split('/')[-1].split(JDi)

hdf5 = name_parts[0] + JDi + '_' + JDf + name_parts[-1]
hdf5 = os.path.join(savepath, hdf5)

print 'Writing:'
print hdf5
uvd.write_uvh5(hdf5, clobber=True)