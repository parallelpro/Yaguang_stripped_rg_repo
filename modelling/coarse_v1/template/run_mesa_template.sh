#!/bin/csh
#PBS -P astero
#PBS -N cs1
#PBS -l mem=31GB
#PBS -l nodes=1:ppn=12
#PBS -l walltime=24:00:00
#PBS -m ea
#PBS -M yali4742@uni.sydney.edu.au
#PBS -V

date
hostname
# module load Anaconda3-5.1.0
source ~/.cshrc


cd "/project/RDS-FSC-astero-RW/low-mass-red-giants/hpc/coarse_v1/template/"
./clean
./mk
python3 "driver_template.py" 1 1


date
exit

###PBS -l nodes=node21:ppn=12
###PBS -l nodes=node43:ppn=12
###PBS -q physics
