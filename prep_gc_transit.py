#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

from __future__ import print_function, division, absolute_import
import numpy as np
import os
import re
from astropy import units
from astropy.time import Time
from pyuvdata import UVData

# define the RA of the galactic center
# RA is 17h45m40.04s in J2000 according to wikipedia
gc_ra = units.Quantity(17 + 45 / 60 + 40.04 / 3600, units.hourangle)

# define width of time centered on the galactic center transit to keep
delta_t = units.Quantity(4 / 60, units.hourangle)
lst_min = gc_ra - delta_t / 2
lst_max = gc_ra + delta_t / 2

# define the polarizations
pols = ['xx', 'yy', 'xy', 'yx']

# loop over miriad files in a directory
jd_dirs = [d for d in sorted(os.listdir(os.getcwd())) if '24575' in d]
for jd in jd_dirs:
    # status update
    print(jd)
    abspath = os.path.abspath(jd)
    basename = os.path.basename(jd)
    os.chdir(jd)

    # make a list of the raw xx data files
    pattern = 'zen\.[0-9]{7}\.[0-9]{5}\.xx\.HH\.uvcRP'
    miriad_files = [f for f in sorted(os.listdir(os.getcwd())) if re.match(pattern, f)]
    for fn in miriad_files:
        uvd = UVData()
        uvd.read_miriad(fn)
        lst_array = units.Quantity(np.unique(uvd.lst_array), units.radian)
        # find out where a sign change happens when differencing lst_array with max and min values;
        # shows where the boundaries are
        dlst_min = lst_array - lst_min
        dlst_max = lst_array - lst_max
        dlst_min_sign = np.sign(dlst_min)
        dlst_max_sign = np.sign(dlst_max)
        signchange_min = ((np.roll(dlst_min_sign, 1) - dlst_min_sign) != 0).astype(int)
        signchange_max = ((np.roll(dlst_max_sign, 1) - dlst_max_sign) != 0).astype(int)
        if np.count_nonzero(signchange_min) > 0:
            min_fn = fn
        if np.count_nonzero(signchange_max) > 0:
            max_fn = fn

    # now go back and read in all the relevant files
    imin = miriad_files.index(min_fn)
    imax = miriad_files.index(max_fn)

    base_list = miriad_files[imin:imax+1]
    fn_list = []
    for fn in base_list:
        obs_list = [re.sub('\.xx\.', '.' + pol + '.', fn) for pol in pols]
        fn_list.extend(obs_list)

    uvd = UVData()
    print("Reading " + ' '.join(fn_list) + '...')
    for i in range(len(fn_list) // 4):
        i1 = 4 * i
        i2 = 4 * (i + 1)
        if i == 0:
            uvd.read_miriad(fn_list[i1:i2])
        else:
            uvd2 = UVData()
            uvd2.read_miriad(fn_list[i1:i2])
            uvd += uvd2
            del(uvd2)
    # only select the relevant times
    indices = np.logical_and(uvd.lst_array > lst_min.to('radian').value,
                             uvd.lst_array < lst_max.to('radian').value)
    print("Selecting times...")
    uvd.select(times=uvd.time_array[indices])

    # write out to a single uvfits file
    outfile = 'gc.' + jd + '.uvcRP.uvfits'
    print("Phasing...")
    uvd.phase_to_time(Time(uvd.time_array[0], format='jd', scale='utc'))
    print("Writing " +  outfile + '...')
    uvd.write_uvfits(outfile, spoof_nonessential=True, clobber=True)

    # clean up
    del(uvd)
    os.chdir('..')
