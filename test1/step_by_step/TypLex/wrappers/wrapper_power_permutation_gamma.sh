#!/bin/bash

qsub -t 0-12   -l walltime=24:00:00,mem=16GB -o ./qsub_out/wrapper_power_permutation_gamma.out -e ./qsub_out/wrapper_power_permutation_gamma.err -v myvar='/home/rf02/rezvan/test1/step_by_step/TypLex/avg_and_stat/Step2_3_MNE_TF_Cluster_Permutation_gamma.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
