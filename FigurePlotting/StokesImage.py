#!/bin/sh

#  Make Stoke Image of Any FITS File using APLPY
#
#
#  Created by Tashalee Billings on 9/12/18.
#

"""
    Change files_from_fits variable to whatever list of fits image files you want.
"""
import os, glob, sys, aplpy
import numpy as np, matplotlib.pyplot as plt

from astropy.io import fits

#-----------------------------------------------
# Create and save .fits plots
#-----------------------------------------------
#***********************************************
files_from_fits = ["imagemodel_gc.2457548.uvc_KBcal.image.fits"] # for real data
#files_from_fits = ["new_gc.2457755.uvC.image.fits"] # for simulated data

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
        
        fig.show_colorscale(cmap=cmap,vmax=vmax,vmin=vmin)
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

    plt.tight_layout()
    fig.savefig('{}.STOKES.png'.format(fitsfile))
    fig.savefig('{}.STOKES.pdf'.format(fitsfile))
    plt.close()
