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

ap = argparse.ArgumentParser(description="Make plots of gain solutions")

ap.add_argument("--fname", type=str, help="name of calfits filt", required=True)
ap.add_argument("--undo_delays", default=False, action="store_true", help="remove delay-induced phase wraps from gains")
ap.add_argument("--kcal", type=str, default=None, help="name of file containing delay solutions")

# define antennas to exclude
bad_ants = [22, 43, 80, 81]

# define output directory
outdir = '/data4/paper/HERA19Golden/kohn18_paper/plots2'

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
    gains = 1j / gains

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
    f1 = plt.figure()
    ax1 = plt.gca()
    f2 = plt.figure()
    ax2 = plt.gca()
    f3 = plt.figure()
    ax3 = plt.gca()
    f4 = plt.figure()
    ax4 = plt.gca()
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
        ax1.plot(freq_array, np.abs(e_gain), linestyle=linestyle, label="{:d}".format(ant))
        ax2.plot(freq_array, np.angle(e_gain), linestyle=linestyle, label="{:d}".format(ant))
        ax3.plot(freq_array, np.abs(n_gain), linestyle=linestyle, label="{:d}".format(ant))
        ax4.plot(freq_array, np.angle(n_gain), linestyle=linestyle, label="{:d}".format(ant))
        counter += 1

    # make plots pretty
    for i, ax in enumerate([ax1, ax2, ax3, ax4]):
        ax.set_xlim((115, 185))
        ax.set_xlabel(r'Frequency [MHz]')
        if i % 2 == 0:
            ax.set_ylim((0, 80))
            ax.set_ylabel(r'Gain amplitude [(Jy/corr unit)$^{1/2}$]')
        else:
            ax.set_ylim((-np.pi, np.pi))
            ax.set_ylabel(r'Gain phase')
        leg = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    # save plots
    output = os.path.join(outdir, 'e_amp.pdf')
    print("Saving {}...".format(output))
    f1.savefig(output, bbox_inches='tight')
    output = os.path.join(outdir, 'e_phs.pdf')
    print("Saving {}...".format(output))
    f2.savefig(output, bbox_inches='tight')
    output = os.path.join(outdir, 'n_amp.pdf')
    print("Saving {}...".format(output))
    f3.savefig(output, bbox_inches='tight')
    output = os.path.join(outdir, 'n_phs.pdf')
    print("Saving {}...".format(output))
    f4.savefig(output, bbox_inches='tight')

    # clean up
    plt.close(f1)
    plt.close(f2)
    plt.close(f3)
    plt.close(f4)

    return

if __name__ == '__main__':
    # parse args
    args = ap.parse_args()

    main(args)
