#!/bin/bash

qsub -t 0-16 -l walltime=24:00:00,mem=16GB -o ./qsub_out/wrapper_filtering.out -e ./qsub_out/wrapper_filtering.err -v myvar='/home/rf02/rezvan/test1/step_by_step/Step1_MNE_Interpolate_Filtering.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
