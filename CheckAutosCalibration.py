#!/bin/sh

#  CheckAutosCalibration.py
#

import numpy as np
import matplotlib.pyplot as plt
from pyuvdata import UVData
import os
from glob import glob
import hera_pspec as hp
from hera_pspec.data import DATA_PATH
import time

sim_path = '/lustre/aoc/projects/hera/plaplant/HERA19Golden/Simulation/'
# Read in the simulated data
sim_autos = UVData()
sim_autos.read_uvh5(sim_path+'zen.2457756.14819.xx.HH.uvCP.autos.uvh5')

# For these purposes, it matters not what the frequency array is, so long as they are all the same
freqs = np.unique(sim_autos.freq_array) # I worry the sims don't quite match the data

# Bring in the beam calibration to K
cosmo = hp.conversions.Cosmo_Conversions()
# List of beamfile to load. This is a healpix map.
beamfile = os.path.join(DATA_PATH, 'NF_HERA_Beams.beamfits')
# intantiate beam and pass cosmology, if not fed, a default Planck cosmology will be assumed
uvb = hp.pspecbeam.PSpecBeamUV(beamfile, cosmo=cosmo)
# find conversion factor from Jy to K.  xx and yy are the same.
Jy_to_K = uvb.Jy_to_mK(freqs,pol='xx')/1000.

# Plot it
plt.plot(freqs,Jy_to_K)
plt.show()

print sim_autos.lst_array.min()
print sim_autos.lst_array.max()
print sim_autos.freq_array.shape

def lstfreqwaterfall(uvd,antval,hdfname):
    plt.figure(figsize=(10,6))
    lst = uvd.lst_array.squeeze()*24./(2.*np.pi)
    freq = uvd.freq_array.squeeze()/1e6
    data = uvd.get_data(bl,force_copy=True)
    bl=(antval,antval)
    flags = ~uvd.get_flags(bl,force_copy=True)
    T = data.real/flags*Jy_to_K
    plt.imshow(T,aspect='auto',
               extent=[freq.min(),freq.max(),lst.max(),lst.min()],vmax=3000,vmin=0)
    plt.ylim([23,10.6])
    plt.colorbar()
    plt.savefig("{}.pdf".format(hdfname))
    return


JDs= ['48','49','50','51','52','53','54','55']
antindex = [9, 10, 20, 22, 31, 43, 53, 64, 65, 72, 80, 81, 88, 89, 96, 97, 104, 105, 112]

for jd in JDs:
    datapath = '/lustre/aoc/projects/hera/plaplant/HERA19Golden/CalibratedData/24575{}/'.format(jd)

    # Now for that juicy real data ...
    realfiles = glob(datapath+'zen.24575{}.*.xx.HH.uvcRPCS.uvh5'.format(jd))

    for ant in antindex:
    
        real_autos = UVData()
        real_autos.read_uvh5(realfiles[0],bls=(ant,ant))
        for realfile in realfiles[1:]:
            print realfile
            tmp = UVData()
            tmp.read_uvh5(realfile,bls=(ant,ant))
            real_autos += tmp
        hdf_name = 'zen.24575{}'.format(jd)+'.xx.HH.auto{}.uvcRPCS.uvh5'.format(ant)
        real_autos.write_uvh5(hdf_name)
        real_autos.read_uvh5(hdf_name)
        
        # Make waterfall plot and save it.
        lstfreqwaterfall(sim_autos,ant,hdfname)
        lstfreqwaterfall(real_autos,ant,hdfname)
