#!/bin/bash

qsub -t 0-14   -l walltime=24:00:00,mem=16GB -o ./qsub_out/wrapper_allinoneconn_permutation.out -e ./qsub_out/wrapper_allinoneconn_permutation.err -v myvar='/home/rf02/rezvan/test1/step_by_step/TypLex/avg_and_stat/Step3_4_MNE_Connectivity_micoh_Permutation.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
