""" This script will save a copy of all given UVP files, averaged over time and baseline.
Reading in the outputs of this script should save time in loading and processing while plotting.
Output files are saved in the same directory with 'avg' appended to their name.
"""

import argparse
import hera_pspec
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument(
    '-F',
    '--files',
    help='Designate the UVP files to be averaged.',
    nargs='*',
    required=True)
args = parser.parse_args()
dfiles = np.array(sorted(args.files))
    
for dfile in dfiles:
    print 'loading:', dfile
    #load files
    uvp = hera_pspec.UVPSpec()
    uvp.read_hdf5(dfile)
    
    print 'averaging:', dfile
    #find baseline groups, average over those groups and time
    blp_grps, lens, angs, tags = hera_pspec.utils.get_blvec_reds(uvp, bl_error_tol=1.0, match_bl_lens=True)
    uvp_avg = uvp.average_spectra(blpair_groups=blp_grps, time_avg=True, inplace=False)
    
    #write out
    uvp_avg.write_hdf5(dfile+'avg', overwrite=True)