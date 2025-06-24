#!/bin/bash
#SBATCH -J LSSTSIMLIBnoHOST
#SBATCH -A m1727
#SBATCH --qos=shared
#SBATCH --constraint=cpu
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8GB
#SBATCH -e ./logs/LSST_SIMLIB_v4.3_%A.err
#SBATCH -o  ./logs/LSST_SIMLIB_v4.3_%A.out


python -u /pscratch/sd/b/bastienc/Soft/OpSimSummaryV2/scripts/make_simlib.py \
'baseline_v4.3_10yrs.db' \
--author 'bastienc' \
-d \
--simlib_coadd 
