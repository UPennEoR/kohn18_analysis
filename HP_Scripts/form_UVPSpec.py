"""This script will take one uvf5 file with a UVData object and turn it into multiple
uvf5 files with UVPspec objects in them that each represent the power spectra of a group
of redundant baselines.  Assumes visibilities are in Jy.
Based on 2_HPWpipe by Austin Fox Fortino.  To be run after select_LSTs_and_concatenate.py."""

import os
import glob
import argparse
import numpy as np
import hera_pspec as hp
from hera_pspec.data import DATA_PATH
from pyuvdata import UVData

parser = argparse.ArgumentParser()
parser.add_argument(
    '-F',
    '--files',
    help='Designate the hdf5 files to be analyzed.',
    nargs='*',
    required=True)
parser.add_argument(
    '-P',
    '--pols',
    help='Designate which polarizations to analyze (e.g.: "pI pQ yx yy xx").',
    nargs='*',
    required=True)
parser.add_argument(
    '-R',
    '--FREQrng',
    help='Designate the frequency range, in channels, that will be analyzed (e.g.: "580 680").',
    type=int,
    nargs=2,
    required=True)
parser.add_argument(
    '-S',
    '--savepath',
    help='Designate the path where the new hdf5 files will be saved. Default is path to data files. A new directory "FREQrng_start_end" will be created.')
args = parser.parse_args()

"""Formatting command line arguments:"""
dfiles = np.array(sorted(args.files))
pols = sorted(args.pols)
FREQrng = args.FREQrng
if args.savepath is None:
    savepath = os.path.dirname(args.files[0])
else:
    savepath = args.savepath
savepath = os.path.join(savepath, 'FREQrng_{}_{}'.format(FREQrng[0], FREQrng[1]))
os.system('mkdir -p {}'.format(savepath))

"""Defining pol constants:"""
STD_POLS = ['xx', 'xy', 'yx', 'yy']
pS_POLS = ['pI', 'pQ', 'pU', 'pV']

"""Loading input files as UVData objects:"""
print "Loading input files"
uvds_std_pols = {std_pol: None for std_pol in STD_POLS}
for dfile in dfiles:
    uvd = UVData()
    uvd.read_uvh5(dfile)
    pol = uvd.get_pols()[0].lower()
    uvds_std_pols[pol] = uvd
    del uvd
    
"""Creating pseudo stokes UVData objects (if requested) and formatting UVdata objects into a dict:"""
uvds = {pol: None for pol in pols}
for pol in pols:
    print "Creating data for pol: " + pol
    
    if pol in STD_POLS:
        uvds[pol] = uvds_std_pols[pol]

    elif pol in pS_POLS:
        if pol == 'pI':
            uvdI = hp.pstokes.construct_pstokes(dset1=uvds_std_pols['xx'], dset2=uvds_std_pols['yy'], pstokes='pI')
            uvds[pol] = uvdI
        if pol == 'pQ':
            uvdQ = hp.pstokes.construct_pstokes(dset1=uvds_std_pols['xx'], dset2=uvds_std_pols['yy'], pstokes='pQ')
            uvds[pol] = uvdQ
        if pol == 'pU':
            uvdU = hp.pstokes.construct_pstokes(dset1=uvds_std_pols['xy'], dset2=uvds_std_pols['yx'], pstokes='pU')
            uvds[pol] = uvdU
        if pol == 'pV':
            uvdV = hp.pstokes.construct_pstokes(dset1=uvds_std_pols['xy'], dset2=uvds_std_pols['yx'], pstokes='pV')
            uvds[pol] = uvdV

"""Making UVPspec objects:"""
for pol, uvd in uvds.items():
    print 'Making UVPSpec object for pol: ' + pol
    # Apply flags
    uvd = hp.flags.construct_factorizable_mask([uvd])[0]
    
    # Intialize a cosmology and a beam
    if pol in STD_POLS:
        beamfile = os.path.join(DATA_PATH, 'HERA_NF_dipole_power.beamfits')
    elif pol in pS_POLS:
        beamfile = os.path.join(DATA_PATH, 'HERA_NF_pstokes_power.beamfits')
    cosmo = hp.conversions.Cosmo_Conversions()
    uvb = hp.pspecbeam.PSpecBeamUV(beamfile, cosmo=cosmo)

    # Convert to cosmological units (mK)
    Jy_to_mK = uvb.Jy_to_mK(np.unique(uvd.freq_array), pol=pol)
    uvd.data_array *= Jy_to_mK[None, None, :, None]

    # Shift data and load datasets
    uvd1 = uvd.select(times=np.unique(uvd.time_array)[:-1:2], inplace=False)
    uvd2 = uvd.select(times=np.unique(uvd.time_array)[1::2], inplace=False)
    ds = hp.PSpecData(dsets=[uvd1, uvd2], wgts=[None, None], beam=uvb)

    # Set visibility units
    ds.dsets[0].vis_units = 'mK'
    ds.dsets[1].vis_units = 'mK'

    # Phase data (What does this do?)
    ds.rephase_to_dset(0)

    """Categorize baselines into physical separation length"""
    xants = uvd.extra_keywords['xants']
    ants = uvd.get_ENU_antpos(pick_data_ants=True)[-1]
    #sort antpairs by their physical separations:
    reds, lens, angs = hp.utils.get_reds(uvd)
    bl_lens = []
    bl_groups = []
    while len(lens) > 0:
        close = np.isclose(lens[0], lens, rtol=0.0, atol=1.0)
        red_lens = np.average(np.take(lens, np.where(close)))
        bl_lens.append(red_lens)
        red_bls = np.concatenate((np.take(reds, np.where(close))[0]))
        bls = []
        for i in range(len(red_bls)):
            ant1, ant2 = red_bls[i]
            if ant1 > ant2:
                ant = ant2
                ant2 = ant1
                ant1 = ant
                del ant
            if (ant1 not in xants) and (ant2 not in xants) and (ant1 in ants) and (ant2 in ants):
                bls.append((ant1, ant2))
        bl_groups.append(bls)
        lens = np.delete(lens, np.where(close))
        reds = np.delete(reds, np.where(close))

    #format norms, bls_reds, baselines, blpairs, and blp_reds to match Austin's conventions
    norms = list(np.round(bl_lens, 1))
    bls_reds = dict(zip(norms, bl_groups))
    baselines = []
    blpairs = []
    for bl_group in bl_groups:
        for bl in bl_group:
            baselines.append(bl)
            blpairs.append((bl, bl))
    blp_reds = {}
    for norm in norms:
        blp_reds[norm] = [(bl, bl) for bl in bls_reds[norm]]

    """Make UVPspec object"""
    uvp = ds.pspec(
        baselines,
        baselines,
        (0, 1),
        pols=[(pol, pol)],
        spw_ranges=[(FREQrng[0], FREQrng[-1])],
        taper="blackman-harris",
        verbose=False)

    """Name and save UVPspec object"""
    end_str = dfiles[0].split('HH')[-1]
    hdf5 = 'zen.{JD}.{JDt0}_{JDtf}.{pol}.HH{end_str}.UVP'.format(
        JD=uvd.extra_keywords['JD'],
        JDt0=uvd.extra_keywords['JDt0'],
        JDtf=uvd.extra_keywords['JDtf'],
        pol=pol,
        end_str=end_str)
    hdf5 = os.path.join(savepath, hdf5)
    print 'Writing: ' + hdf5
    uvp.write_hdf5(hdf5, overwrite=True)
    
"""Saving additional metadata:"""
print "Saving additional metadata"
metadata = os.path.join(savepath, 'metadata.npz')
np.savez(
    metadata,
    uvd_extra_keywords=uvd.extra_keywords,
    bls_reds=bls_reds,
    blp_reds=blp_reds,
    baselines=baselines,
    blpairs=blpairs,
    norms=norms,
    antpos={ant: pos for ant, pos in zip(uvd.get_ENU_antpos()[1], uvd.get_ENU_antpos()[0])},
    FREQrng=FREQrng)