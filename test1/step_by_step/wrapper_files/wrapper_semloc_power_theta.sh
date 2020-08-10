#!/bin/bash

qsub -t 0-16  -l walltime=24:00:00,mem=16GB -o ./qsub_out/wrapper_semloc_power_theta.out -e ./qsub_out/wrapper_semloc_power_theta.err -v myvar='/home/rf02/rezvan/test1/step_by_step/Step6_1_MNE_SourceBand_power_theta.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
