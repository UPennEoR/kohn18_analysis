#!/bin/bash
# Copyright (c) 2018 UPennEoR
# Licensed under the BSD 2-clause License

# This script is calibrate_data.sh. It is intended to calibrate all of the
# HERA19 Golden Set data with the same bandpass calibration and absolute scale
# reference spectrum.

# activate conda environment
source activate hera

# define directory where data live
RAWDATA_DIR=/data4/paper/HERA19Golden/RawData
CALDATA_DIR=/data4/paper/HERA19Golden/CalibratedData

# define calibration solution files
# we use just bandpass solutions and absolute ratio scaling
BCAL=/data4/paper/HERA19Golden/CalibratedData/calibration_solutions/imagemodel_gc.2457548.uvcRP.MS.B.cal.h5
ACAL=/data4/paper/HERA19Golden/CalibratedData/calibration_solutions/abscal_spec.h5

# define location of python script for generating calfits files
CALFITS_SCRIPT=/data4/paper/HERA19Golden/kohn18_analysis/hdf5_to_calfits.py

# loop over data directories
for JD in 2457548 2457549 2457550 2457551 2457552 2457553 2457554 2457555; do
    # move to raw data directory
    cd ${RAWDATA_DIR}/${JD}

    # make a list of all xx files
    xx_files=`ls *.xx.HH.uvcRP.uvh5`
    for xx_fn in ${xx_files}; do
	# make a calfits file
        calfits_fn=`echo ${xx_fn} | sed s/uvh5/calfits/ | sed s/\.xx\././`
        echo python ${CALFITS_SCRIPT} --bcal ${BCAL} --acal ${ACAL} --overwrite --fname ${calfits_fn} --uv_file ${xx_fn}
	python ${CALFITS_SCRIPT} --bcal ${BCAL} --acal ${ACAL} --overwrite --fname ${calfits_fn} --uv_file ${xx_fn}

	# apply the calibration to each polarization
	xy_fn=`echo ${xx_fn} | sed s/\.xx\./.xy./`
	yx_fn=`echo ${xx_fn} | sed s/\.xx\./.yx./`
	yy_fn=`echo ${xx_fn} | sed s/\.xx\./.yy./`
	xx_fn_out=`echo ${xx_fn} | sed s/uvcRP/uvcRPC/`
	xy_fn_out=`echo ${xy_fn} | sed s/uvcRP/uvcRPC/`
	yx_fn_out=`echo ${yx_fn} | sed s/uvcRP/uvcRPC/`
	yy_fn_out=`echo ${yy_fn} | sed s/uvcRP/uvcRPC/`
	echo apply_cal.py --new_cal ${calfits_fn} --filetype_in uvh5 --filetype_out uvh5 --gain_convention divide --clobber ${xx_fn} ${xx_fn_out}
	apply_cal.py --new_cal ${calfits_fn} --filetype_in uvh5 --filetype_out uvh5 --gain_convention divide --clobber ${xx_fn} ${xx_fn_out}
	echo apply_cal.py --new_cal ${calfits_fn} --filetype_in uvh5 --filetype_out uvh5 --gain_convention divide --clobber ${xy_fn} ${xy_fn_out}
	apply_cal.py --new_cal ${calfits_fn} --filetype_in uvh5 --filetype_out uvh5 --gain_convention divide --clobber ${xy_fn} ${xy_fn_out}
	echo apply_cal.py --new_cal ${calfits_fn} --filetype_in uvh5 --filetype_out uvh5 --gain_convention divide --clobber ${yx_fn} ${yx_fn_out}
	apply_cal.py --new_cal ${calfits_fn} --filetype_in uvh5 --filetype_out uvh5 --gain_convention divide --clobber ${yx_fn} ${yx_fn_out}
	echo apply_cal.py --new_cal ${calfits_fn} --filetype_in uvh5 --filetype_out uvh5 --gain_convention divide --clobber ${yy_fn} ${yy_fn_out}
	apply_cal.py --new_cal ${calfits_fn} --filetype_in uvh5 --filetype_out uvh5 --gain_convention divide --clobber ${yy_fn} ${yy_fn_out}

	# move output files to calibrated directory
	mv ${xx_fn_out} ${CALDATA_DIR}/${JD}
	mv ${xy_fn_out} ${CALDATA_DIR}/${JD}
	mv ${yx_fn_out} ${CALDATA_DIR}/${JD}
	mv ${yy_fn_out} ${CALDATA_DIR}/${JD}
    done
done
