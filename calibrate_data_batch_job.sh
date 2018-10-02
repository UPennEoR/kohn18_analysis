#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -l h_vmem=12G
#$ -j y
#$ -N caldata
#$ -o log
#$ -l all
#$ -m be
#$ -M plaplant@sas.upenn.edu

# call script to do calibration
./calibrate_data.sh 2>&1 caldata.log
