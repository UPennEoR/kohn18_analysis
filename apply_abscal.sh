#!/bin/bash
# Copyright (c) 2018 UPennEoR
# Licensed under the BSD 2-clause License

# This script is apply_abscal.sh. It is used to apply a calibration solution
# saved in a .calfits file to a uvdata set. This script expects positional
# arguments to specify the calibration file, input uv dataset, and output
# uv dataset.
#
# Usage:
#   apply_abscal.sh <calfits_file> <input_uv> <output_uv>
#
# Notes:
#   The input and output filetypes are inferred from the extensions.
#   If the extension is "uvfits" or "uvh5", the filetype is set accordingly.
#   All other extensions are assumed to be miriad files.

# infer input filetype from name of input
input_ext="${2##*.}"
if [ input_ext == "uvfits" ]; then
    filetype_in="uvfits"
else
    filetype_in="miriad"
fi

# infer output filetype from name of output
output_ext="${3##*.}"
if [ output_ext == "uvfits" ]; then
    filetype_out="uvfits"
elif [ output_ext == "uvh5" ]; then
    filetype_out="uvh5"
else
    filetype_out="miriad"
fi

apply_cal.py --new_cal "${1}" --filetype_in "${filetype_in}" --filetype_out "${filetype_out}" --clobber "${2}" "${3}"
