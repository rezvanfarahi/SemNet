#!/bin/bash

qsub -t 0-17 -l walltime=24:00:00,mem=32GB -o ./qsub_out/wrapper_evoked_permutation.out -e ./qsub_out/wrapper_evoked_permutation.err -v myvar='/home/rf02/rezvan/semnet/avg_and_stat/Step3_1_MNE_Evoked_Cluster_Permutation.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
