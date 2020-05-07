#!/bin/bash

qsub -t 0-4,7-8,10-13,15-24 -l walltime=24:00:00,mem=32GB -o ./qsub_out/wrapper_invop.out -e ./qsub_out/wrapper_invop.err -v myvar='/home/rf02/rezvan/semnet/newfiltscripts/Step5_MNE_Noise_Covariance_AND_Inverse_Operator.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
