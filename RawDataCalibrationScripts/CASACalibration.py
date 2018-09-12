#!/bin/sh

#  SaulPaperCASA.py
#  
#
#  Created by Tashalee Billings on 9/6/18.
#  

"""
   The REAL DATA has the correct antenna positions added and the uvw added to the array (so you do not need to run add_uvw.py). Simply combine the 4 polarized miriad files and [xx,yy,xy,yx] and write to uvfits.
   
   This script will go calibrate and image the data.
   
   1) KB calibrate to a single flat-spectrum 1 Jy pt src at the GC,
   2) do an interactive CLEAN to get a better GC model,
   3) repeat the KB calibration using the CLEAN model, and
   4) write down the gain solutions for the KB CLEAN model,
   5) dirty image the calibrated data and spectrum image
   
"""

# Location of the data on folio2
#/data4/paper/HERA19Golden/RawData/24575*/*uvcR
#/data4/paper/HERA19Golden/RawData/{JD}/gc.{JD}.uvcRP.uvfits 4-min long GC centered files
# zen.2457548.46619.HH.uvcR.MS

import numpy as np
import time
import shutil
import argparse
import os, glob, sys


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

def imagemodel(namems, mi): #NAMEofIMAGE.model
    ft(namems,model=mi,usescratch = True) # multiple flat point spectrum model.
    return

def flatspecmodel(namems): # single flat point spectrum model.
    cl.addcomponent(flux=1.0, fluxunit='Jy', shape='point',
                    dir='J2000 17h45m40.0409s -29d0m28.118s')
    
    if os.path.exists("GC.cl"):
        shutil.rmtree("GC.cl")
    cl.rename("GC.cl")
    cl.close()
    
    ft(namems, complist="GC.cl", usescratch=True)
    return

def specmodel(namems): # Frequency sprectrum model
    direction = 'J2000 17h45m40.0409s -29d0m28.118s'
    ref_freq = '408MHz'

    cl.addcomponent(label="Galactic Center", flux=3709, fluxunit='Jy',
                    dir=direction, freq=ref_freq, shape='point',
                    spectrumtype='spectral index',index=-0.5)
    
    if os.path.exists("GC_spec.cl"):
        shutil.rmtree("GC_spec.cl")
    cl.rename("GC_spec.cl")
    cl.close()
                    
    ft(namems, complist="GC_spec.cl", usescratch=True)
    return

def delay(namems):
    calitable = os.path.basename(namems) + ".K.cal"
    gaincal(namems,caltable=calitable,gaintype = 'K' , solint='inf',refant='10')
    print("Created "+glob.glob(calitable)[0])
    applycal(namems,gaintable=[calitable])

    return

def Gphase(namems):
    calitable = os.path.basename(namems) + ".Gphase.cal"
    gaincal(namems,caltable=calitable,gaintype = 'G', solint='int', refant='10', calmode='p')
    print("Created "+glob.glob(calitable)[0])
    applycal(namems,gaintable=[calitable])
    
    return


def Gampl(namems):
    calitable = os.path.basename(namems) + ".Gamp.cal"
    gaincal(namems,caltable=calitable,gaintype = 'G', solint='int', refant='10', calmode='a')
    print("Created "+glob.glob(calitable)[0])
    applycal(namems,gaintable=[calitable])
    
    return

def bandpasscal(namems):
    calitable = os.path.basename(namems) +".B.cal"
    bandpass(namems,caltable=calitable,solint='inf',combine='scan',refant='10')
    print("Created "+glob.glob(calitable)[0])
    applycal(namems,gaintable=[calitable])
    
    return

def clean1(namems):
    clean(namems,namems[:-5]+'_no_cal_niter0',niter=0, weighting = 'briggs',robust =0, imsize =[512 ,512] ,pbcor=False, cell=['500 arcsec'] ,mode='mfs',nterms=1, spw='0:150~900', stokes='IQUV',interactive=False, npercycle=5, threshold='0.1mJy/beam')

    return

parser = argparse.ArgumentParser(description='Perform Delay and Bandpass Calibration. Run with casa as: casa -c SaulPaperCASA.py --<args>=')

parser.add_argument("-c", type=str, help="Calibrationg for different permutations.")
#parser.add_argument("--refant", type=str, default=None, help="The reference Antenna used to calibrate.")

#   command_line = "casa --nologger --nocrashreport --nogui --agg -c ../../calibration_tests_in_CASA.py --refant=10 --nocal='False' --calmode_G='p'"

args = parser.parse_args()

if __name__ == '__main__': #This only runs is we run directly but not if you try to import it.
    # Convert from uvfits to MS
    uvfits = glob.glob('*.uvfits')
    for uvfit in uvfits:
        msfile = 'flatspecmodel_'+uvfit.strip('uvfits') + 'MS'
        importuvfits(vis=msfile,fitsfile=uvfit)
        print(glob.glob("flatspecmodel_*.MS")[0]+" created.")
    
    msname_ = glob.glob("flatspecmodel_gc.*.MS")

    for msname in msname_:
        
        print("The measurement set you are working with is ... "+msname)
        time.sleep(2)
        # Single Flat Spectrum Point Source Calibration
        flag(msname)
        flatspecmodel(msname)

        print("Beginning The Calibration Process.")
        # Start calibrtion Process
        delay(msname) # Delay Calibration
        bandpasscal(msname) # Bandpass Calibration

        print("Interactive CLEAN")
        time.sleep(2)
        # Do interactive CLEAN to create image.model file that will be used as a model to calibrate to.
        clean(msname,msname[:-5]+'_KGcal_model',niter=0, weighting = 'briggs',robust =0, imsize =[512 ,512] ,pbcor=False, cell=['500 arcsec'] ,mode='mfs',nterms=1, spw='0:150~900', stokes='IQUV',interactive=True, npercycle=5, threshold='0.1mJy/beam')

        modelname = glob.glob(msname[:-5]+'_KGcal_model.model')[0]
        
        # Multiple Flat Spectrum Point Source Calibration
        
        for uvfit in uvfits:
            msfile = "imagemodel_"+uvfit.strip('uvfits') + "MS"
            importuvfits(vis=msfile,fitsfile=uvfit)
            print(glob.glob("imagemodel_*.MS")[0]+" created.")

        msnamei = glob.glob("imagemodel_*.MS")[0] # rename MS

        print("The measurement set you are working with is ... "+msnamei)
        time.sleep(2)
        flag(msnamei)
        imagemodel(msnamei, mi=modelname)

        print("Beginning The Calibration Process.")
        # Start calibrtion Process
        delay(msnamei) # Delay Calibration
        bandpasscal(msnamei) # Bandpass Calibration

        cparamB = glob.glob("*.B.cal")
        fparamK = glob.glob("*.K.cal")
        
        for bc in cparamB:
            tb.open(bc)
            bp = tb.getcol('CPARAM') # (2,1024,113) was orginally (1024,113)
            bp_ants = tb.getcol('ANTENNA1') # (113,)
            bp_flags = tb.getcol('FLAG') # (2,1024,113) was orginally (1024,113)
            tb.close()
            # load spectral window data
            tb.open(bc+"/SPECTRAL_WINDOW")
            bp_freqs = tb.getcol('CHAN_FREQ') #(1024, 1)
            tb.close()
            # write to file
            np.savez("{}.npz".format(bc), bp=bp, bp_ants=bp_ants, bp_flags=bp_flags, bp_freqs=bp_freqs)
            print("Created "+ bc+'.npz')

        for kc in fparamK:
            tb.open(kc)
            delays = tb.getcol('FPARAM') # (2, 1, 113) was orginally (113,)
            delay_ants = tb.getcol('ANTENNA1') # (113,)
            delay_flags = tb.getcol('FLAG') # (2, 1, 113) was orginally (113,)
            tb.close()
            np.savez("{}.npz".format(kc), delay_ants=delay_ants, delays=delays, delay_flags=delay_flags)
            print("Created "+ kc+'.npz')

        print("CLEAN")
        # Dirty Image of Calibrated Data
        clean(msnamei,msnamei[:-5]+'_KBcal',niter=0, weighting = 'briggs',robust =0, imsize =[512 ,512] ,pbcor=False, cell=['500 arcsec'] ,mode='mfs',nterms=1, spw='0:150~900', stokes='IQUV',interactive=False, npercycle=5, threshold='0.1mJy/beam')
        exportfits(msnamei[:-5]+'_KBcal.image',msnamei[:-5]+'_KBcal.image.fits')

        # Make Spectrum Image
        fitsname = 'mult_spec_chan_'+msnamei[:-3]+'.image.fits'
        imagename = 'mult_spec_chan_'+msnamei[:-3]+'.image'
    
        clean(msnamei,imagename[:-6],niter=0,weighting = 'briggs',robust =0,
              imsize =[512 ,512],cell=['500 arcsec']
              ,mode='channel',nterms =1,spw='0:150~900',nchan=750, start=150, width=1
              ,stokes='IQUV', interactive=False, npercycle=5, threshold='0.1mJy/beam')
        exportfits(imagename,fitsname)

