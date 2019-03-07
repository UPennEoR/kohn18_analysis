"""Plots 4 UVPSpec objects side by side.  Saves both 1D and 2D plots of the data."""

import os
import argparse
import numpy as np
import matplotlib
matplotlib.use('agg')
matplotlib.rc('text', usetex=True)
matplotlib.rc('font', family='serif')
import matplotlib.pyplot as plt
import hera_pspec
import pyuvdata

parser = argparse.ArgumentParser()
parser.add_argument(
    '-D',
    '--data_files',
    help='Designate the real data UVP files to be plotted.  Expects 4 files, one for each pI, pQ, pU, pV.',
    nargs='*',
    required=True)
parser.add_argument(
    '-S',
    '--sim_files',
    help='Designate the sim data UVP files to be plotted.  Expects 4 files, one for each pI, pQ, pU, pV.',
    nargs='*',
    required=True)
parser.add_argument(
    '-P',
    '--path',
    help='Designate the path where the new hdf5 files will be saved. Default is path to data files.')
parser.add_argument(
    '-N',
    '--name',
    help='Designate the filename for the plots, without the extension.',
    required=True)
parser.add_argument(
    '-E',
    '--extension',
    help='Designate the extension to save the plots as.  Default is "pdf".',
    default='pdf')
parser.add_argument(
    '-V',
    '--vminmax',
    help='Designate vmin and vmax in log10 for the 2D plots.',
    required=True,
    nargs=2,
    type=float)
parser.add_argument(
    '-Q',
    '--quiet',
    help='Toggle off print statements',
    action='store_false')
args = parser.parse_args()

#parse arguments
dfiles = np.array(sorted(args.data_files))
sfiles = np.array(sorted(args.sim_files))
name = args.name
ext = args.extension
vmin, vmax = args.vminmax
verbose = args.quiet
if args.path is None:
    savepath = os.path.dirname(args.files[0])
else:
    savepath = args.path

#load UVPSpec objects
duvps = []
suvps = []
pols = ['pI','pQ','pU','pV']
if verbose:
    print 'reading in observed data files'
for dfile in dfiles:
    uvp = hera_pspec.UVPSpec()
    uvp.read_hdf5(dfile)
    duvps.append(uvp)
if verbose:
    print 'reading in simulated data files'
for sfile in sfiles:
    uvp = hera_pspec.UVPSpec()
    uvp.read_hdf5(sfile)
    suvps.append(uvp)

#plot 1D plot
def plot_1D_row(fig, axes, uvps, pols, labeltop=False, labelbot=False, scalexax=False):
    for ax, uvp, pol in zip(axes, uvps, pols):
        #form the subplot
        blp_grps, lens, angs, tags = hera_pspec.utils.get_blvec_reds(uvp, bl_error_tol=1.0, match_bl_lens=True)
        f = hera_pspec.plot.delay_spectrum(uvp, blp_grps, 0, pol, average_blpairs=True, average_times=True, legend=False, ax=ax, delay=False)
        if pol != 'pI':
            ax.set_ylabel('')
        if pol == 'pV':
            ax.legend([str(np.round(l,1))+' m' for l in lens])
            
        axT = ax.twiny()
        delays = uvp.get_dlys(0) * 10**9
        axT.set_xlim([delays[0],delays[-1]])
        axT.set_xlabel(r"$\tau$ $[{\rm ns}]$", fontsize=16, labelpad=10)

        scale = .465
        axT.set_xlim([lim*scale for lim in axT.get_xlim()])
        if pol=='pV' and scalexax:
            ax.set_xlim([lim*scale for lim in ax.get_xlim()])
        
        if not labeltop:
            axT.set_xlabel('')
            axT.set_xticks([])
        else:
            ax.set_title(pol,pad=45)  
        if not labelbot:
            ax.set_xlabel('')

        horizons = uvp.get_blpair_seps() / hera_pspec.conversions.units.c * 1e9
        horizons = np.array([horizons[-1]])
        avg_z = uvp.cosmo.f2z(np.mean(uvp.freq_array[uvp.spw_to_freq_indices(0)]))
        horizons *= uvp.cosmo.tau_to_kpara(avg_z) / 1e9
        for horizon in horizons:
            ax.axvline(horizon, linestyle='--', color='k', linewidth=.7)
            ax.axvline(-horizon, linestyle='--', color='k', linewidth=.7)

        ax.grid()

#plot 2D plot
def plot_2D_row(fig, axes, uvps, pols, vmin, vmax, cbax, labeltop=False, labelbot=False, scalexax=False, plotcbar=False):
    for ax, uvp, pol in zip(axes, uvps, pols):
        blp_grps, lens, angs, tags = hera_pspec.utils.get_blvec_reds(uvp, bl_error_tol=1.0, match_bl_lens=True)
        colorbar = pol=='pV' and plotcbar
        if pol=='pI':
            ax.invert_yaxis()
        f = hera_pspec.plot.delay_wedge(uvp, 0, pol, center_line=True, horizon_lines=True, colorbar=colorbar, cbax=cbax, lw=1, flip_yax=False, ax=ax, vmin=vmin, vmax=vmax, delay=False)
        if pol != 'pI':
            ax.set_ylabel('')
        if pol=='pV':
            axR = ax.twinx()
            axR.invert_yaxis()
            axR.set_ylim(ax.get_ylim())
            axR.set_yticks(ax.get_yticks())
            axR.set_yticklabels(np.round(lens,1), fontsize=11)
            axR.set_ylabel(r"$|\vec{b}|$ $[{\rm m}]$", fontsize=16)
        axT = ax.twiny()
        delays = uvp.get_dlys(0) * 10**9
        axT.set_xlim([delays[0],delays[-1]])
        axT.set_xlabel(r"$\tau$ $[{\rm ns}]$", fontsize=16, labelpad=10)

        if not labeltop:
            axT.set_xlabel('')
            axT.set_xticks([])
        else:
            ax.set_title(pol,pad=45)  
        if not labelbot:
            ax.set_xlabel('')
        
        scale = .48
        axT.set_xlim([lim*scale for lim in axT.get_xlim()])
        if pol=='pV' and scalexax:
            ax.set_xlim([lim*scale for lim in ax.get_xlim()])

#ready subplots and colorbar axes
fig, axes = plt.subplots(3, 4, figsize=(15, 18), sharex=True, sharey='row')
cbax = fig.add_axes([.88,0,.1,1])
cbax.set_frame_on(False)
cbax.set_xticks([])
cbax.set_yticks([])
plt.subplots_adjust(wspace=0, hspace=0)

#plot each row
if verbose:
    print 'plotting row 1 of 3'
plot_2D_row(fig=fig, axes=axes[0], uvps=suvps, vmin=vmin, vmax=vmax, cbax=cbax, pols=pols, labeltop=True)
if verbose:
    print 'plotting row 2 of 3'
plot_2D_row(fig=fig, axes=axes[1], uvps=duvps, vmin=vmin, vmax=vmax, cbax=cbax, pols=pols, plotcbar=True)
if verbose:
    print 'plotting row 3 of 3'
plot_1D_row(fig=fig, axes=axes[2], uvps=duvps, pols=pols, labelbot=True, scalexax=True)
#since sharex=True for all graphs, only scalexax=True once at the end, otherwise it gets scaled multiple times

#save the fig
fname = os.path.join(savepath, '{name}.{ext}'.format(name=name, ext=ext))
if verbose:
    print 'saving figure to: ' + fname
fig.savefig(fname, bbox_inches='tight')