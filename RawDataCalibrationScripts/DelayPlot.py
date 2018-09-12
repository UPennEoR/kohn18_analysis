#!/bin/sh

#  SaulPaperDelayPlot.py
#  
#
#  Created by Tashalee Billings on 9/11/18.
#  
import numpy as np, matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.colors as mcolors
import os, glob, sys, aplpy
import time
import h5py

from astropy import units as u
from astropy.io import fits
from matplotlib import gridspec

Knpzlist = glob.glob("*.K.cal.npz") #(2,1,113)
antindex = [9, 10, 20, 22, 31, 43, 53, 64, 65, 72, 80, 81, 88, 89, 96, 97, 104, 105, 112] #position of ant list?

#***********************************************
#-----------------------------------------------
# Create and save .NPZ plots
#-----------------------------------------------
#***********************************************

for kc in Knpzlist:
    k=np.load(kc)["delays"]
    badant = np.where(k[0,0,:]==0.0)[0]
    k[0,0,badant]='nan'
    k[1,0,badant]='nan'
    
    plt.figure(figsize=(14,8))
    
    plt.plot(k[0,0,:],'.')
    plt.title("Delay Values: East Pol")
    plt.xlabel('Antenna Number')
    plt.ylabel('Amp [ns]')
    plt.grid(which='major', linestyle='-', linewidth='0.5', color='red')
    plt.grid(which='minor', linestyle='-', linewidth='0.5', color='black')
    
    plt.savefig("DELAY_amp_east_{}.png".format(kc))
    plt.close()
    
    plt.figure(figsize=(14,8))
    plt.plot(k[1,0,:],'.')
    plt.title("Delay Values: North Pol")
    plt.xlabel('Antenna Number')
    plt.ylabel('Amp [ns]')
    plt.grid(which='major', linestyle='-', linewidth='0.5', color='red')
    plt.grid(which='minor', linestyle='-', linewidth='0.5', color='black')
    
    plt.savefig("DELAY_amp_north_{}.png".format(kc))
    plt.close()


