#!/bin/bash
#PBS -N cal_H19
#PBS -q hera
#PBS -l nodes=1:ppn=1
#PBS -l vmem=8G
#PBS -l walltime=48:00:00
#PBS -j oe
#PBS -o log
#PBS -V
#PBS -m be
#PBS -M plaplant@sas.upenn.edu

# call script to do calibration
cd $PBS_O_WORKDIR
./calibrate_data_nrao.sh 2>&1 caldata.log
