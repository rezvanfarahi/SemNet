#!/bin/bash

qsub -t 0-15   -l walltime=24:00:00,mem=16GB -o ./qsub_out/wrapper_tempconnectivity_permutation.out -e ./qsub_out/wrapper_tempconnectivity_permutation.err -v myvar='/home/rf02/rezvan/test1/step_by_step/averaging_and_stat/Step3_2_MNE_Connectivity_Temporal_Permutation.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
