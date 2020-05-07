#!/bin/bash

qsub -t 0-20 -l walltime=24:00:00,mem=32GB -o ./qsub_out/wrapper_invop.out -e ./qsub_out/wrapper_invop.err -v myvar='/home/rf02/rezvan/semnet/semnet4semloc/Step5_MNE_Noise_Covariance_AND_Inverse_Operator.py'  /home/rf02/rezvan/semnet/semnet4semloc/wrappers/wrapper_batch.sh
