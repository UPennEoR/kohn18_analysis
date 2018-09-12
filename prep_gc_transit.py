#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
# Copyright (c) 2018 UPennEoR
# Licensed under the 2-clause BSD License

from __future__ import print_function, division, absolute_import
import numpy as np
import os
import re
from astropy import units
from astropy.time import Time
from pyuvdata import UVData

# define the directories where data live
real_data = False
if real_data:
    # actual raw data
    data_dirs = ['/data4/paper/HERA19Golden/RawData/2457548',
                 '/data4/paper/HERA19Golden/RawData/2457548',
                 '/data4/paper/HERA19Golden/RawData/2457548',
                 '/data4/paper/HERA19Golden/RawData/2457548',
                 '/data4/paper/HERA19Golden/RawData/2457548',
                 '/data4/paper/HERA19Golden/RawData/2457548',
                 '/data4/paper/HERA19Golden/RawData/2457548',
                 '/data4/paper/HERA19Golden/RawData/2457548']
    uv_ext = 'uvcRP'
else:
    data_dirs = ['/data4/paper/HERA19Golden/Simulation']
    uv_ext = 'uvC'

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
for data_dir in data_dirs:
    # status update
    print(data_dir)
    os.chdir(data_dir)

    # make a list of the raw xx data files
    pattern = 'zen\.[0-9]{{7}}\.[0-9]{{5}}\.xx\.HH\.{uv_ext}'.format(uv_ext=uv_ext)
    miriad_files = [f for f in sorted(os.listdir(os.getcwd())) if re.match(pattern, f)]

    # extract integer JD
    pattern = 'zen\.([0-9]{7})\.'
    jd = re.match(pattern, miriad_files[0]).group(1)

    # initialize filenames for max and min lst (defined above)
    min_fn = None
    max_fn = None
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

    # get the indices corresponding to max and min lst
    # handle edge-cases of transition happening on a boundary; assume that we only need neighboring files
    if min_fn is None and max_fn is None:
        raise ValueError("no sign change for max and min LST ranges; either a "
                         "pathological case of both sign changes falling on file "
                         "boundaries, or no sign changes present in the data.")
    if min_fn is not None:
        imin = miriad_files.index(min_fn)
    else:
        imin = miriad_files.index(max_fn) - 1
    if max_fn is not None:
        imax = miriad_files.index(max_fn)
    else:
        imax = miriad_files.index(min_fn) + 1

    # flesh out list to all 4 polarizations
    base_list = miriad_files[imin:imax+1]
    fn_list = []
    for fn in base_list:
        obs_list = [re.sub('\.xx\.', '.' + pol + '.', fn) for pol in pols]
        fn_list.extend(obs_list)

    # get around weird concatenation bug by combining polarizations first, then
    # adjacent times
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
    outfile = 'gc.' + jd + '.{}.uvfits'.format(uv_ext)
    print("Phasing...")
    uvd.phase_to_time(Time(uvd.time_array[0], format='jd', scale='utc'))
    print("Writing " +  outfile + '...')
    uvd.write_uvfits(outfile, spoof_nonessential=True)

    # clean up
    del(uvd)
