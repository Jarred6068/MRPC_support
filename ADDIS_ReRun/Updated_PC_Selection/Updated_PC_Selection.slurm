#!/bin/bash
#SBATCH -o stdout.out
#SBATCH -e stderr.out

module load R/4.0.2

cd /mnt/ceph/jarredk/AddisReRunFiles/Updated_PC_Selection/

R CMD BATCH Updated_PC_Selection.R 
