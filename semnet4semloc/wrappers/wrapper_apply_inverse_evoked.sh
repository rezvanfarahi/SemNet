#!/bin/bash

qsub -t 0-19 -l walltime=24:00:00,mem=32GB -o ./qsub_out/wrapper_apply_inverse_evoked.out -e ./qsub_out/wrapper_apply_inverse_evoked.err -v myvar='/home/rf02/rezvan/semnet/semnet4semloc/Step6_1_MNE_Applyinverse_evoked.py'  /home/rf02/rezvan/semnet/semnet4semloc/wrappers/wrapper_batch.sh
