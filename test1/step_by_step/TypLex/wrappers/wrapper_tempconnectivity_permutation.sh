#!/bin/bash

qsub -t 0-15   -l walltime=24:00:00,mem=16GB -o ./qsub_out/wrapper_specconnectivity_permutation.out -e ./qsub_out/wrapper_specconnectivity_permutation.err -v myvar='/home/rf02/rezvan/test1/step_by_step/TypLex/avg_and_stat/Step4_3_MNE_WBConnectivity_Temporal_Permutation.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
