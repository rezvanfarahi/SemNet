#!/bin/bash

qsub -t 0-24 -l walltime=24:00:00,mem=16GB -o ./qsub_out/wrapper_filtering.out -e ./qsub_out/wrapper_filtering.err -v myvar='/home/rf02/rezvan/semnet/semnet4semloc/Step2_MNE_Interpolate_Filtering.py'  /home/rf02/rezvan/semnet/semnet4semloc/wrappers/wrapper_batch.sh
