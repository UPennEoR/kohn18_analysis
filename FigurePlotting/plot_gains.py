#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2018 UPennEoR
# Licensed under the 2-clause BSD License

from __future__ import print_function, division, absolute_import
from pyuvdata import UVCal
import numpy as np
import argparse
import os
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

ap = argparse.ArgumentParser(description="Make plots of gain solutions")

ap.add_argument("--fname", type=str, help="name of calfits filt", required=True)
ap.add_argument("--undo_delays", default=False, action="store_true", help="remove delay-induced phase wraps from gains")
ap.add_argument("--kcal", type=str, default=None, help="name of file containing delay solutions")

# define antennas to exclude
bad_ants = [22, 43, 80, 81]

# define channels to exclude
bad_chans = np.zeros(1024, dtype=np.bool)
bad_chans[0:66] = True
bad_chans[930:] = True
bad_chans[377:388] = True
bad_chans[850:855] = True
bad_chans[831] = True
bad_chans[769] = True
bad_chans[511] = True
bad_chans[913] = True
bad_chans[695:705] = True
bad_chans[770:780] = True

# define output directory
outdir = '/data4/paper/HERA19Golden/kohn18_paper/Plots'

def main(args):
    # matplotlib options
    matplotlib.rc('text', usetex=True)
    matplotlib.rc('font', family='serif')

    # read in calfits file
    uvc = UVCal()
    uvc.read_calfits(args.fname)

    # get some metadata
    # hera_cal put duplicates along antenna axis, so only take every other antenna
    ant_array = uvc.ant_array[::2]
    freq_array = uvc.freq_array[0, :] * 1e6  # convert to Hz; accidentally saved as MHz
    jones_array = uvc.jones_array

    # get just the first time, since all are the same
    # also take reciprocal so we get gains with units of (Jy/corr unit)**(1/2)
    gains = uvc.gain_array[::2, 0, :, 0, :].squeeze()
    gains = (1 + 0j) / gains

    if args.undo_delays:
        # remove phase wraps from delays
        if args.kcal is None:
            raise ValueError("when --undo_gains is passed in, must also specify the gain file with --kcal argument")
        # read in as hdf5 file
        delays = np.zeros_like(ant_array, dtype=np.float)
        with h5py.File(args.kcal, 'r') as f:
            delay_ants = f["/Data/delay_ants"].value.tolist()
            delays = f["/Data/delays"].value.T
            delay_flags = f["/Data/delay_flags"].value.T
        # convert from ns -> s
        delays *= 1e-9

        # reorder antennas
        delays = np.array(map(lambda a: delays[delay_ants.index(a), :] if a in delay_ants else 0., ant_array))

        # convert delays into complex gains; freqs has units of Hz
        delay_gains = np.exp(-2j * np.pi * np.einsum('a,bc->bac', freq_array - freq_array.min(), delays))

        # undo delays from gains
        gains *= delay_gains

    # convert from Hz -> MHz
    freq_array /= 1e6

    # make plots of gain amplitude/phase for e/n pols
    f = plt.figure()
    gs = gridspec.GridSpec(2, 2, hspace=0., wspace=0.1)
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1], sharey=ax1)
    ax3 = plt.subplot(gs[2], sharex=ax1)
    ax4 = plt.subplot(gs[3], sharex=ax2, sharey=ax3)
    counter = 0
    for i, ant in enumerate(ant_array):
        if ant in bad_ants:
            continue
        e_gain = gains[i, :, 0]
        n_gain = gains[i, :, 1]
        if counter >= 10:
            linestyle = '--'
        else:
            linestyle = '-'

        # mask arrays based on bad channels
        abs_e = np.ma.masked_where(bad_chans, np.abs(e_gain))
        abs_n = np.ma.masked_where(bad_chans, np.abs(n_gain))
        phs_e = np.ma.masked_where(bad_chans, np.angle(e_gain))
        phs_n = np.ma.masked_where(bad_chans, np.angle(n_gain))
        ax1.plot(freq_array, abs_e, linestyle=linestyle, label="{:d}".format(ant))
        ax2.plot(freq_array, abs_n, linestyle=linestyle, label="{:d}".format(ant))
        ax3.plot(freq_array, phs_e, linestyle=linestyle, label="{:d}".format(ant))
        ax4.plot(freq_array, phs_n, linestyle=linestyle, label="{:d}".format(ant))
        counter += 1

    # make plots pretty
    ax1.set_xlim((115, 185))
    ax2.set_xlim((115, 185))
    ax1.set_ylim((0, 80))
    ax3.set_ylim((-np.pi, np.pi))
    ax1.set_ylabel(r'Gain amplitude [(Jy/corr unit)$^{1/2}$]')
    ax3.set_ylabel(r'Gain phase')
    ax3.set_xlabel(r'Frequency [MHz]')
    ax4.set_xlabel(r'Frequency [MHz]')
    leg = ax4.legend(loc='center left', bbox_to_anchor=(1, 1))
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    ax1.text(0.8, 0.95, 'E pol', transform=ax1.transAxes, bbox=props,
             verticalalignment='top')
    ax2.text(1.9, 0.95, 'N pol', transform=ax1.transAxes, bbox=props,
             verticalalignment='top')
    ax2.get_yaxis().set_visible(False)
    ax4.get_yaxis().set_visible(False)

    # add shaded regions showing low and high band
    for ax in [ax1, ax2, ax3, ax4]:
        ax.axvspan(120, 130, alpha=0.3, color='k')
        ax.axvspan(157, 167, alpha=0.3, color='k')

    # save plots
    output = os.path.join(outdir, 'gains.pdf')
    print("Saving {}...".format(output))
    f.savefig(output, bbox_inches='tight')

    # clean up
    plt.close(f)

    return

if __name__ == '__main__':
    # parse args
    args = ap.parse_args()

    main(args)
