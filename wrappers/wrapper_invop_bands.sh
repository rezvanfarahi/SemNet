#!/bin/bash

qsub -t 0-19 -l walltime=24:00:00,mem=32GB -o ./qsub_out/wrapper_invop_bands.out -e ./qsub_out/wrapper_invop_bands.err -v myvar='/home/rf02/rezvan/semnet/newfiltscripts/Step5_MNE_Noise_Covariance_AND_Inverse_Operator_bands.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
