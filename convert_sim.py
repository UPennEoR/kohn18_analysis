from pyuvdata import UVData
import numpy as np
from astropy import units as u
from astropy import constants as c
from glob import glob

sim_path = '/data4/paper/SimulatedData/HERA19_GSM2008_NicCST2016_full_day/true_visibilities/'
simfiles = glob(sim_path+'*.uv')

golden_path = '/data4/paper/HERA19Golden/Simulation/'

for simfile in simfiles:
    filename = simfile.split('/')[-1]
    uvd = UVData()
    print 'Reading '+simfile
    uvd.read_miriad(simfile)
    # Fix conjugation error
    uvd.data_array *= -1.
    # Fix calibration to Jy
    # Look, Ma! No loops!
    wvl = (c.c/(uvd.freq_array*u.Hz)).to(u.m).value
    Ksr2Jy = np.reshape(np.outer(np.ones(uvd.Nblts),1.e26*2.*c.k_B.value/np.power(wvl,2).squeeze()),uvd.data_array.shape)
    uvd.data_array *= Ksr2Jy
    outfile = golden_path+filename+'C'
    print 'Writing '+ outfile
    uvd.write_miriad(outfile)
