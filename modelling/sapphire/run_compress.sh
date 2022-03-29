#!/bin/csh
#PBS -P astero
#PBS -N compress
#PBS -l select=1:ncpus=12:mem=60GB
#PBS -l walltime=2:00:00
#PBS -m ea
#PBS -M yali4742@uni.sydney.edu.au
#PBS -V

date
hostname
#module load Anaconda3-5.1.0
cd "/project/RDS-FSC-astero-RW/numax-sc-metallicity/hpc/coarse_v1/"
source "~/.cshrc"
python3 "compress.py" 12


date
exit

###PBS -l nodes=node21:ppn=12
###PBS -l nodes=node43:ppn=12
###PBS -q physics

