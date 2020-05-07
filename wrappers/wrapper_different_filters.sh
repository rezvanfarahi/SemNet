#!/bin/bash

qsub -t 0-18 -l walltime=24:00:00,mem=16GB -o ./qsub_out/wrapper_different_filters.out -e ./qsub_out/wrapper_different_filters.err -v myvar='/home/rf02/rezvan/semnet/Step2_MNE_Choose_highpass_Filter.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
