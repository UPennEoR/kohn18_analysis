#!/bin/sh

#  SaulPaperSpectScalePython.py
#  
#
#  Created by Tashalee Billings on 9/12/18.
#  

"""
   Make Stokes Image for real and simulated data visibilities. Make stokes image of the Primary Beam. Makes a plot overlaying the Spectrum of the Real and Sim data you get from making channel by channel images.
"""

import numpy as np, matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.colors as mcolors
import os, glob, sys, aplpy

from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename

parser = argparse.ArgumentParser(description='Creates Unscaled Flux Density Spectrum. Run with Python as: python flux_desntiy_spectrum.py --<args>=')
parser.add_argument("--fits_image", type=str, default=None, help="Fits image file name. eg 'zen.2457755.89889.HH.uv.TrueVis.image.fits' ")
parser.add_argument("--scaled_fits_image", type=str, default=None, help="Fits image file name. eg 'scaled_zen.2457755.89889.HH.uv.TrueVis.image.fits' ")
parser.add_argument("--sim_fits_image", type=str, default=None, help="Fits image file name. eg 'zen.2457755.89889.HH.uv.TrueVis.image.fits' ")
parser.add_argument("--fits_multspec", type=str, default=None, help="Mult Spec Fits image file name. eg 'mult_spec_chan_zen.2457755.89889.HH.uv.TrueVis.image.fits' ")
parser.add_argument("--sim_fits_multspec", type=str, default=None, help="Mult Spec Fits image file name. eg 'mult_spec_chan_zen.2457755.89889.HH.uv.TrueVis.image.fits' ")
parser.add_argument("--ra_pixel", type=int, default=None, help="Open CASA image file, pick a channel and identify the RA pixel number associated with the largest Flux Density of the GC.")
parser.add_argument("--sim_ra_pixel", type=int, default=None, help="Open CASA image file, pick a channel and identify the RA pixel number associated with the largest Flux Density of the GC.")
parser.add_argument("--dec_pixel", type=int, default=None, help="Open CASA image file, pick a channel and identify the DEC pixel number associated with the largest Flux Density of the GC.")
parser.add_argument("--sim_dec_pixel", type=int, default=None, help="Open CASA image file, pick a channel and identify the DEC pixel number associated with the largest Flux Density of the GC.")

args = parser.parse_args()

if __name__ == '__main__': #This only runs is we run directly but not if you try to import it.
    #ra, dec = 269,251
    #sim_ra,sim_dec = 268,259
    spect_realdata = fits.open(args.fits_multspec) #(4, 1024, 512, 512)
    spect_simdata = fits.open(args.sim_fits_multspec) #(4, 1024, 512, 512)
    
    realdata = fits.open(args.fits_image) # (4,1, 512, 512)
    simdata = fits.open(args.sim_fits_image) # (4,1, 512, 512)
    
    ra, dec = args.ra_pixel,args.dec_pixel
    sim_ra,sim_dec = args.sim_ra_pixel, args.sim_dec_pixel

    #Find the spectrum need to scale the visibilities
    spect_factor = spect_simdata[0].data[0,:,sim_ra,sim_dec]/spect_realdata[0].data[0,:,ra,dec] # (512,)
    nan_index = np.argwhere(np.isnan(spect_factor))[:,0].tolist() # look for any nans
    spect_factor[nan_index]=0 # If they exist set them to zero
    np.savez("SpectrumScale.npz", SpectrumScale =spect_factor)# Save them as a separate file

    #factor = np.mean(np.ma.masked_equal(spect_factor,0)) #(1,)
    
    realdata[0].data *= spect_factor # (4, 1, 512, 512)
    realdata.writeto(args.scaled_fits_image)

    #-----------------------------------------------
    # Create and save .fits plots
    #-----------------------------------------------
    #***********************************************
    files_from_fits = [args.scaled_fits_image,args.sim_fits_image]

    for fitsfile in files_from_fits:
        
        f = plt.figure(figsize=(10,8))
        for pol in np.arange(4):
            fig = aplpy.FITSFigure(fitsfile,dimensions=[0,1],slices=[0,pol],figure=f,subplot=(2,2,pol+1))
            if pol == 0:
                vmax=2220
                vmin=0
                cmap="viridis"
            if pol == 1:
                vmax=20
                vmin=-20
                cmap="RdYlGn"
            if pol == 2:
                vmax=20
                vmin=-20
                cmap="RdYlGn"
            if pol == 3:
                vmax=20
                vmin=-20
                cmap="RdYlGn"
            
            fig.show_colorscale(cmap=cmap)#,vmax=vmax,vmin=vmin)#,stretch='arcsczdxcinh')
            fig.add_grid()
            fig.grid.set_color('black')
            fig.grid.set_xspacing(15)
            fig.grid.set_yspacing(15)
            fig.grid.show()
            fig.axis_labels.set_font(size='small')
            fig.tick_labels.set_font(size='small')
            fig.tick_labels.set_xformat('hh')
            fig.tick_labels.set_yformat('dd')
            fig.add_colorbar()
            fig.colorbar.set_font(size='small')
        
        plt.suptitle('{} STOKE Visibilities'.format(fitsfile))
        plt.tight_layout()
        fig.savefig('{}.STOKES.png'.format(fitsfile))
        fig.savefig('{}.STOKES.pdf'.format(fitsfile))
        plt.close()

# Image of the Simulated Primary Beam
    fitsfile = glob.glob("*.psf.fits")[0]

    f = plt.figure(figsize=(10,8))
    fig = aplpy.FITSFigure(fitsfile,dimensions=[0,1],slices=[0,0],figure=f)
    fig.show_colorscale(cmap='viridis',vmax=1.0,vmin=-0.2)#,stretch='arcsczdxcinh')
    fig.add_grid()
    fig.grid.set_color('black')
    fig.grid.set_xspacing(15)
    fig.grid.set_yspacing(15)
    fig.grid.show()
    fig.axis_labels.set_font(size='small')
    fig.tick_labels.set_font(size='small')
    fig.tick_labels.set_xformat('hh')
    fig.tick_labels.set_yformat('dd')
    fig.add_colorbar()
    fig.colorbar.set_font(size='small')

    #plt.suptitle('{} Primary Beam'.format(fitsfile))
    plt.tight_layout()
    fig.savefig('{}.PrimaryBeam.png'.format(fitsfile))
    fig.savefig('{}.PrimaryBeam.pdf'.format(fitsfile))
#   plt.show()

# Image of the Simulation and Real Data Spectrum
    nimage = np.shape(spect_simdata[0].data)[1]
    freq = np.linspace(115,188,nimage) #(start freq, end freq, num of points between)

    plt.figure(figsize=(15,8))
        
    plt.plot(freq,np.ma.masked_equal(_datasims[0,:,sim_ra,sim_dec],0),'r.', label='Simulation')
    plt.plot(freq,np.ma.masked_equal(_datareal[0,:,ra,dec],0),'b.', label='Real Data')
    plt.xlabel("Freq [Mz]")
    plt.ylabel("Flux Density [Jy/Beam] ")

    #plt.title('{} Flux Density Spectrum'.format(args.sim_fits_multspec))
    plt.legend(loc='upper right')
    plt.savefig('SimAndReal_FluxDensitySpectrum.png')
    plt.savefig('SimAndReal_FluxDensitySpectrum.pdf')
    #plt.show()
    plt.close()
