#!/bin/sh

#  SaulPaperSimData.py
#  
#
#  Created by Tashalee Billings on 9/12/18.
#

"""
   This script only works if it live in the same directory as the data objects.
"""
import os, glob, sys

#Flag the same RFI and bad antennas as Real Data
def flag(namems): #You might have to update this. Flag MS or calibration Tables.
    flagdata(namems, flagbackup=True, mode='manual',antenna="23" ) # 22 for HERA19
    flagdata(namems, flagbackup=True, mode='manual',antenna="44" ) # 43 for HERA19
    flagdata(namems, flagbackup=True, mode='manual',antenna="81" ) # 80 for HERA19
    flagdata(namems, flagbackup=True, mode='manual',antenna="82" ) # 81 for HERA19
    
    flagdata(namems, flagbackup=True, mode='manual', spw="0:0~65" )#channels 0-65 of spw 0
    flagdata(namems, flagbackup=True, mode='manual', spw ="0:930~1023")
    flagdata(namems, flagbackup=True, mode='manual', spw="0:377~387" )
    flagdata(namems, flagbackup=True, mode='manual', spw="0:850~854" )
    flagdata(namems, flagbackup=True, mode='manual', spw = "0:831" )
    flagdata(namems, flagbackup=True, mode='manual', spw = "0:769" )
    flagdata(namems, flagbackup=True, mode='manual', spw = "0:511" )
    flagdata(namems, flagbackup=True, mode='manual', spw = "0:913" )
    
    flagdata(namems, autocorr = True )
    return

#convert uvfits to MS
uvfits = glob.glob('*.uvfits')
for uvfit in uvfits:
    msfile = uvfit.strip('uvfits') + 'MS'
    importuvfits(vis=msfile,fitsfile=uvfit)
    print(glob.glob("*.MS")[0]+" created.")

msname_ = glob.glob("*.MS")

for msname in msname_:
    flag(msname)

    # Make dirty image of Simulated data and export to fits
    imgname = msname[:-5]+'_TrueVis_NoDeconv'
    fitsimg = msname[:-5]+'_TrueVis_NoDeconv.image.fits'

    clean(msname,imgname,niter=0, weighting = 'briggs',robust =0, imsize =[512 ,512] ,pbcor=False,
          cell=['500 arcsec'] ,phasecenter = 'J2000 17h45m40.0409s -29d0m28.118s',
          mode='mfs',nterms=1, spw='0:150~900', stokes='IQUV',
          interactive=False, npercycle=5, threshold='0.1mJy/beam')
    exportfits(imgname+'.image',fitsimg)
    
    # Make Interactive CLEAN image of Simulated data and export to fits
    imgnameC = msname[:-5]+'_TrueVis'
    fitsimgC = msname[:-5]+'_TrueVis.image.fits'
    
    clean(msname,imgnameC,niter=0, weighting = 'briggs',robust =0, imsize =[512 ,512] ,pbcor=False,
          cell=['500 arcsec'],phasecenter = 'J2000 17h45m40.0409s -29d0m28.118s',
          mode='mfs',nterms=1, spw='0:150~900', stokes='IQUV',interactive=True, npercycle=5, threshold='0.1mJy/beam')
    exportfits(imgnameC+'.image',fitsimgC)

    # Make Spectrum Image and export to fits
    fitsname = 'mult_spec_1024_chan_'+msname[:-3]+'.image.fits'
    imagename = 'mult_spec_1024_chan_'+msname[:-3]+'.image'
    
    clean(msname,imagename[:-6],niter=0,weighting = 'briggs',robust =0
          ,imsize =[512 ,512],cell=['500 arcsec'],phasecenter = 'J2000 17h45m40.0409s -29d0m28.118s'
          ,mode='channel',nterms =1,spw='0',nchan=1024, start=0, width=1
          ,stokes='IQUV', interactive=False, npercycle=5, threshold='0.1mJy/beam')
    exportfits(imagename,fitsname)



