# -*- coding: utf-8 -*-
# Copyright (c) 2018 UPennEoR
# Licensed under the 2-clause BSD License

from __future__ import print_function, division, absolute_import
from pyuvdata import UVData
import numpy as np
from astropy import units
from astropy import constants as const
from glob import glob

sim_path = '/data4/paper/SimulatedData/HERA19_GSM2008_NicCST2016_full_day/true_visibilities/'
simfiles = glob(sim_path + '*.uv')

golden_path = '/data4/paper/HERA19Golden/Simulation/'

for simfile in simfiles:
    filename = simfile.split('/')[-1]
    uvd = UVData()
    print('Reading ' + simfile)
    uvd.read_miriad(simfile)
    # Fix conjugation error
    uvd.data_array = np.conj(uvd.data_array)
    # Fix calibration to Jy
    # Look, Ma! No loops!
    wvl = (const.c / (uvd.freq_array*units.Hz)).to(units.m).value
    # Factor of 2 should be checked.  This might disagree with hera_pspec convention for
    # forming Stokes parameters.  I think it is independent of the beam, because it is applied as 
    # a conversion from a sky model in Kelvin to one in Jy/sr, and can be pulled outside the 
    # visibility integral.  
    Ksr2Jy = np.reshape(np.outer(np.ones(uvd.Nblts), 1.e26 * 2. * const.k_B.value
                                 / np.power(wvl,2).squeeze()), uvd.data_array.shape)
    uvd.data_array *= Ksr2Jy
    uvoutfile = golden_path + filename+'C'
    print('Writing ' + uvoutfile)
    uvd.write_miriad(uvoutfile,clobber=True)
    uvh5outfile = golden_path+filename+'C.uvh5'
    print('Writing '+uvh5outfile)
    uvd.write_uvh5(uvh5outfile,clobber=True)
