#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2018 UPennEoR
# Licensed under the 2-clause BSD License

import aipy as a, numpy as np, scipy.constants as c
import optparse, sys, os
import pyuvdata

o = optparse.OptionParser()
o.set_usage('add_uvws.py [options] *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o,cal=True)
opts,args = o.parse_args(sys.argv[1:])

print 'opts.cal',opts.cal

#uv_latlonalt = pyuvdata.UVData()
#uv_latlonalt.read_miriad('zen.2458098.60274.xx.HH.uv')
#lat, lon, alt = uv_latlonalt.telescope_location_lat_lon_alt
lat, lon, alt = (-0.5361918109651213, 0.37399448506783717, 1073.000000008382)

for uvfile in args:
    uvofile = uvfile+'U'
    print uvfile,'->',uvofile
    if os.path.exists(uvofile):
        print uvofile, 'exists, skipping.'
        continue
    
    uvi = a.miriad.UV(uvfile)
    uvo = a.miriad.UV(uvofile, status='new')
    
    #aa = a.phs.ArrayLocation(('-30:43:17.5','21:25:41.9'))
    #get the positions in the rotated ECEF coordinates
    aa = a.cal.get_aa(opts.cal,np.array([0.15]))
    ENU_antpos = aa.antpos_ideal.T
    ECEF_antpos = pyuvdata.ECEF_from_ENU(ENU_antpos, lat, lon, alt)
    rotECEF_antpos = pyuvdata.rotECEF_from_ECEF(ECEF_antpos.T, lon) * 10**9 / c.c
    nints = 0
    curtime = None
    
    #format it so unused ants are at 0, 0, 0, as pyuvdata expects
    for i in range(len(ENU_antpos[0])):
        if (ENU_antpos.T[i]==[-1.,-1.,-1.]).all():
            rotECEF_antpos[i] = [0.,0.,0.]
    
    #make sure that there are the appropriate number of antennae in the position file
    nants = uvi['nants']
    cal_nants = rotECEF_antpos.shape[0]
    if nants > cal_nants:
        for i in range(nants-cal_nants):
            rotECEF_antpos = np.append(rotECEF_antpos, [[0.,0.,0.]], axis=0)
    elif nants < cal_nants:
        rotECEF_antpos = rotECEF_antpos[:nants]
    
    #flatten coords into a single list, as miriad expects them to be saved
    rotECEF_antpos = rotECEF_antpos.T.flatten()
    
    def mfunc(uv, preamble, data, flags):
        global curtime
        global nints
        #print curtime, nints
        uvw, t, (i,j) = preamble
        #print 'uvw',uvw
        del uvw
        uvw = aa.get_baseline(i,j,'z')
        #print 'uvw from aa'
        #print uvw
        preamble = (uvw, t, (i,j))
        return preamble, data, flags

    uvo.init_from_uv(uvi,override={'antpos':rotECEF_antpos})
    uvo.pipe(uvi, mfunc=mfunc, raw=True, append2hist=' '.join(sys.argv)+'\n')
