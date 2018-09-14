#!/bin/sh

#  SaulPaperSpectScalePython.py
#
#
#  Created by Tashalee Billings on 9/12/18.
#

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
    # SIM 268,259 or # REAL 269,251
    spect_realdata = fits.open(args.fits_multspec) #(4,750, 512, 512)  from chan 150-900
    spect_simdata = fits.open(args.sim_fits_multspec) #(4,750, 512, 512) from chan 150-900
    realdata = fits.open(args.fits_image) # (4,1, 512, 512)
    simdata = fits.open(args.sim_fits_image) # (4,1, 512, 512)
    
    ra, dec = args.ra_pixel,args.dec_pixel
    sim_ra,sim_dec = args.sim_ra_pixel, args.sim_dec_pixel
    
    spect_factor = spect_simdata[0].data[0,:,sim_ra,sim_dec]/spect_realdata[0].data[0,:,ra,dec] #(750,) Spectrum of GC for Stokes I
    nans = np.argwhere(np.isnan(spect_factor)).tolist()
    spect_factor[nans] = 0
    meanval = np.mean(np.ma.masked_equal(spect_factor,0)) #(1,) 5129.22761707989
    
    realdata[0].data *=*meanval # (4,1, 512, 512)
    realdata.writeto(args.scaled_fits_image)
    
    #-----------------------------------------------
    # Create and save .fits plots
    #-----------------------------------------------
    #***********************************************
    files_from_fits = [args.scaled_fits_image,args.sim_fits_image]
    
    for fitsfile in files_from_fits:
        
        f = plt.figure(figsize=(15,8))
        for pol in np.arange(4):
            fig = aplpy.FITSFigure(fitsfile,dimensions=[0,1],slices=[0,pol],figure=f,subplot=(2,2,pol+1))
            if pol == 0:
                vmax=1.7
                vmin=-0.2
                cmap="viridis"
            if pol == 1:
                vmax=0.3
                vmin=-.2
                cmap="RdYlGn"
            if pol == 2:
                vmax=0.1
                vmin=-0.03
                cmap="RdYlGn"
            if pol == 3:
                vmax=0.003
                vmin=-0.007
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
        fig.savefig('{}.STOKES.png'.format(fitsfile))
        plt.close()
    
    # Image of the Simulated Primary Beam
    fitsfile = "new_gc.2457755.u_TrueVis.psf.fits"
    
    f = plt.figure(figsize=(15,8))
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
    fig.savefig('{}.PrimaryBeam.png'.format(fitsfile))
    plt.show()

