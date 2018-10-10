"""Plots 4 UVPSpec objects side by side.  Saves both 1D and 2D plots of the data."""

import os
import argparse
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import hera_pspec
import pyuvdata

parser = argparse.ArgumentParser()
parser.add_argument(
    '-F',
    '--files',
    help='Designate the UVP files to be plotted.  Expects 4 files, one for each pI, pQ, pU, pV.',
    nargs='*',
    required=True)
parser.add_argument(
    '-W',
    '--wedge',
    help='Toggle to turn 2D wedge plotting off',
    action='store_false',
    default=True)
parser.add_argument(
    '-O',
    '--one',
    help='Toggle to turn 1D wedge plotting off',
    action='store_false',
    default=True)
parser.add_argument(
    '-S',
    '--savepath',
    help='Designate the path where the new hdf5 files will be saved. Default is path to data files.')
parser.add_argument(
    '-N',
    '--name',
    help='Designate the filename for the plots.  "1D" or "2D" will be appended.  Default is an empty string.',
    default='')
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
args = parser.parse_args()

#parse arguments
dfiles = np.array(sorted(args.files))
dfiles_basename = [os.path.basename(dfile) for dfile in dfiles]
wedge = args.wedge
one = args.one
name = args.name
if len(name) > 0:
    name = name + '_'
ext = args.extension
vmin, vmax = args.vminmax
if args.savepath is None:
    savepath = os.path.dirname(args.files[0])
else:
    savepath = args.savepath

#load UVPSpec objects
uvps = []
pols = ['pI','pQ','pU','pV']
for dfile in dfiles:
    uvp = hera_pspec.UVPSpec()
    uvp.read_hdf5(dfile)
    uvps.append(uvp)

#plot 1D plot
if one:
    #ready the subplots
    fig, axes = plt.subplots(1, 4, figsize=(15, 6), sharey=True, sharex=True)
    plt.subplots_adjust(wspace=0, hspace=0.1)
    
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
        axT.set_xlabel(r"$\tau$ $[{\rm ns}]$", fontsize=16)

        scale = .465
        axT.set_xlim([lim*scale for lim in axT.get_xlim()])
        if pol=='pV':
            ax.set_xlim([lim*scale for lim in ax.get_xlim()])

        horizons = uvp.get_blpair_seps() / hera_pspec.conversions.units.c * 1e9
        horizons = np.array([horizons[-1]])
        avg_z = uvp.cosmo.f2z(np.mean(uvp.freq_array[uvp.spw_to_freq_indices(0)]))
        horizons *= uvp.cosmo.tau_to_kpara(avg_z) / 1e9
        for horizon in horizons:
            ax.axvline(horizon, linestyle='--', color='k', linewidth=.7)
            ax.axvline(-horizon, linestyle='--', color='k', linewidth=.7)

        ax.set_title(pol,pad=45)
        ax.grid()

    plt.savefig(os.path.join(savepath, '{name}1D.{ext}'.format(
        name=name,
        ext=ext)),
        bbox_inches='tight')

#plot 2D plot
if wedge:
    #ready subplots and colorbar axes
    fig, axes = plt.subplots(1, 4, figsize=(15, 6), sharex=True, sharey=True)
    cbax = fig.add_axes([.88,0,.1,1])
    cbax.set_frame_on(False)
    cbax.set_xticks([])
    cbax.set_yticks([])
    plt.subplots_adjust(wspace=0, hspace=0.1)

    for ax, uvp, pol in zip(axes, uvps, pols):
        blp_grps, lens, angs, tags = hera_pspec.utils.get_blvec_reds(uvp, bl_error_tol=1.0, match_bl_lens=True)
        colorbar=False
        if pol=='pI':
            ax.invert_yaxis()
        if pol=='pV':
            colorbar=True
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
        axT.set_xlabel(r"$\tau$ $[{\rm ns}]$", fontsize=16)

        scale = .48
        axT.set_xlim([lim*scale for lim in axT.get_xlim()])
        if pol=='pV':
            ax.set_xlim([lim*scale for lim in ax.get_xlim()])

        ax.set_title(pol,pad=45)
        
    plt.savefig(os.path.join(savepath, '{name}2D.{ext}'.format(
        name=name,
        ext=ext)),
        bbox_inches='tight')
