#!/bin/bash

qsub -t 0-19 -l walltime=24:00:00,mem=32GB -o ./qsub_out/wrapper_apply_inverse_evoked.out -e ./qsub_out/wrapper_apply_inverse_evoked.err -v myvar='/home/rf02/rezvan/semnet/newfiltscripts/Step6_1_MNE_Applyinverse_evoked_bands.py'  /home/rf02/rezvan/test1/wrapper_batch.sh