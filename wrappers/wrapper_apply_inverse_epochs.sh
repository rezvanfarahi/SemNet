#!/bin/bash

qsub -t 0-18 -l walltime=24:00:00,mem=64GB -o ./qsub_out/wrapper_apply_inverse_epochs.out -e ./qsub_out/wrapper_apply_inverse_epochs.err -v myvar='/home/rf02/rezvan/semnet/Step8_MNE_Applyinverse_epochs.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
