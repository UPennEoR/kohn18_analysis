#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2018 UPennEoR
# Licensed under the 2-clause BSD License

from __future__ import print_function, division, absolute_import
from pyuvdata import UVData, UVCal
import numpy as np
import argparse
import os
import h5py
import hera_cal

ap = argparse.ArgumentParser(description="Convert CASA gain solutions in hdf5 files to .calfits files")

ap.add_argument("--fname", type=str, help="name of calfits file", required=True)
ap.add_argument("--uv_file", type=str, help="name of uvfile to apply calibration to", required=True)
ap.add_argument("--kcal", type=str, default=None, help="name of file containing delay solutions")
ap.add_argument("--bcal", type=str, default=None, help="name of file containing bandpass solutions")
ap.add_argument("--acal", type=str, default=None, help="name of file containing absolute scale reference spectrum")
ap.add_argument("--overwrite", default=False, action="store_true", help="overwrite output file if it exists")
ap.add_argument("--multiply_gains", default=False, action="store_true", help="change gain convention from divide to multiply")
ap.add_argument("--smooth_ratio", default=False, action="store_true", help="smooth sim/data ratio with a low-order polynomial")


def main(ap):
    # read in UVData object
    uv_file = args.uv_file
    uvd = UVData()
    uvd.read(uv_file)

    # get metadata
    antpos, ants = uvd.get_ENU_antpos(center=True, pick_data_ants=True)
    ants = list(ants)
    freqs = np.unique(uvd.freq_array[0, :])
    times = np.unique(uvd.time_array)
    Nants = len(ants)
    Nfreqs = len(freqs)
    Ntimes = len(times)
    # assume we have linear polarizations x and y for gain solutions, in that order
    jones_array = ['x', 'y']
    Njones = len(jones_array)

    # start with unity gains
    gains = np.ones((Nants, Nfreqs, 1, Njones), dtype=np.complex)
    flags = np.zeros((Nants, Nfreqs, 1, Njones), dtype=np.bool)

    # make sure the user specified at least one gain solution
    if args.kcal is None and args.bcal is None and args.acal is None:
        raise ValueError("at least one calibration type must be specified")

    # read in K-gain file
    kfile = args.kcal
    if kfile is not None:
        delays = np.zeros_like(ants, dtype=np.float)
        with h5py.File(kfile, 'r') as f:
            delay_ants = f["/Data/delay_ants"].value.tolist()
            delays = f["/Data/delays"].value.T
            delay_flags = f["/Data/delay_flags"].value.T
        # convert from ns -> s
        delays *= 1e-9

        # reorder antennas
        delays = np.array(map(lambda a: delays[delay_ants.index(a), :] if a in delay_ants else 0., ants))
        delay_flags = np.array(map(lambda a: delay_flags[delay_ants.index(a), :] if a in delay_ants else True, ants))

        # convert delays into complex gains; freqs has units of Hz
        delay_gains = np.exp(-2j * np.pi * np.einsum('a,bc->bac', freqs - freqs.min(), delays))[:, :, np.newaxis, :]
        delay_gain_flags = delay_flags[:, np.newaxis, :]
        delay_gain_flags = np.repeat(delay_gain_flags, Nfreqs, axis=1)[:, :, np.newaxis, :]

        assert(np.allclose(delay_gains[0, :, 0, 0], np.exp(-2j * np.pi * (freqs - freqs.min()) * delays[0, 0])))
        assert(np.allclose(delay_gains[0, :, 0, 1], np.exp(-2j * np.pi * (freqs - freqs.min()) * delays[0, 1])))
        assert(np.allclose(delay_gains[-1, :, 0, 0], np.exp(-2j * np.pi * (freqs - freqs.min()) * delays[-1, 0])))
        assert(np.allclose(delay_gains[-1, :, 0, 1], np.exp(-2j * np.pi * (freqs - freqs.min()) * delays[-1, 1])))

        # make sure we have the right shape
        right_shape = (Nants, Nfreqs, 1, Njones)
        if delay_gains.shape != right_shape:
            raise ValueError("bandpass gains are not the right shape; expecting {}, got {}".format(right_shape, delay_gains.shape))

        # multiply into gains
        gains *= delay_gains

        # add flags
        flags += delay_gain_flags

    # read in B-gain file
    bfile = args.bcal
    if bfile is not None:
        with h5py.File(bfile, 'r') as f:
            # casa saves bandpasses as (Njones, Nfreqs, Nants)
            # we want to transpose them to (Nants, Nfreqs, Njones)
            bp_gains = np.swapaxes(f["/Data/bp"].value, 0, 2)
            bp_ants = f["/Data/bp_ants"].value.tolist()
            bp_freqs = f["/Data/bp_freqs"].value
            bp_flags = np.swapaxes(f["/Data/bp_flags"].value, 0, 2)

        # get number of frequencies
        bp_Nfreqs = len(bp_freqs)

        # reorder antennas
        bp_gains = np.array([bp_gains[bp_ants.index(a), :, :] if a in bp_ants
                             else np.ones((bp_Nfreqs, Njones), dtype=np.complex) for a in ants])
        bp_flags = np.array([bp_flags[bp_ants.index(a), :, :] if a in bp_ants
                             else np.ones((bp_Nfreqs, Njones), dtype=np.bool) for a in ants])

        # get gains and flags into the right shape
        bp_gains = bp_gains[:, :, np.newaxis, :]
        bp_flags = bp_flags[:, :, np.newaxis, :]

        # make sure we have the right shape
        right_shape = (Nants, Nfreqs, 1, Njones)
        if bp_gains.shape != right_shape:
            raise ValueError("bandpass gains are not the right shape; expecting {}, got {}".format(right_shape, bp_gains.shape))

        # multipy into gains
        gains *= bp_gains.conj()

        # add flags
        flags += bp_flags

    # read in overall amplitude spectrum
    afile = args.acal
    if afile is not None:
        filename, ext = os.path.splitext(afile)
        if ext == '.npz':
            f = np.load(afile)
            amp = f["SpectrumScale"]
        elif ext == '.h5' or ext == '.hdf5':
            # assume spectrum is stored in hdf5 format
            with h5py.File(afile, 'r') as f:
                amp = f["/Data/spectrum_scale"].value
        else:
            raise ValueError("unrecognized filetype for abscale spectrum")

        # optionally smooth the data
        if args.smooth_ratio:
            # define frequencies -- assume 1024 frequency channels between 100 and 200 MHz
            freqs = np.linspace(100., 200., num=1024, endpoint=False)
            # define cutoff channels for where to compute polynomial fit
            freq_low = 150
            freq_hi = 900

            # generate 3rd order polynomial
            p = np.poly1d(np.polyfit(freqs[freq_low:freq_hi], amp[freq_low, freq_hi], 3))
            amp = p(freqs)

        # turn it into the right shape
        amp = np.stack((amp,) * Njones).T
        amp = np.stack((amp,) * Nants)
        right_shape = (Nants, Nfreqs, Njones)
        if amp.shape != (Nants, Nfreqs, Njones):
            raise ValueError("amp shape is {}; was expecting {}".format(amp.shape, right_shape))
        assert(np.allclose(amp[0, :, 0], amp[1, :, 0]))
        assert(np.allclose(amp[0, :, 0], amp[0, :, 1]))

        # add an extra dummy time axis
        amp = amp[:, :, np.newaxis, :]

        # multiply into the gains; we take the square root so that g_i * g_j^* gives the original
        # also, some frequencies have values of inf, so we need to handle those separately
        sqrt_amp = np.where(amp == np.inf, np.inf, np.sqrt(amp))
        gains *= sqrt_amp

        # adjust flags to account for inf values
        amp_flags = np.where(amp == np.inf, True, False).astype(np.bool)
        flags += amp_flags

    # make the gains and flags the right number of time samples
    gains = np.repeat(gains, Ntimes, axis=2)
    flags = np.repeat(flags, Ntimes, axis=2)

    # make into calfits
    fname = args.fname
    overwrite = args.overwrite
    if args.multiply_gains:
        gain_convention = 'multiply'
    else:
        gain_convention = 'divide'
    gain_dict = {}
    flag_dict = {}
    for i, ant in enumerate(ants):
        for j, pol in enumerate(jones_array):
            gain_dict[(ant, pol)] = gains[i, :, :, j].T
            flag_dict[(ant, pol)] = flags[i, :, :, j].T
    uvc = hera_cal.io.write_cal(fname, gain_dict, freqs, times, flags=flag_dict,
                                overwrite=overwrite, gain_convention=gain_convention)

    return


if __name__ == '__main__':
    # parse args
    args = ap.parse_args()

    main(ap)
