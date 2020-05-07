#!/bin/bash

qsub -t 0-10 -l walltime=24:00:00,mem=32GB -o ./qsub_out/wrapper_signedevoked_permutation_vsbaseline.out -e ./qsub_out/wrapper_signedevoked_permutation_vsbaseline.err -v myvar='/home/rf02/rezvan/semnet/avg_and_stat/Step3_1_MNE_signedEvoked_and_Cluster_Permutation_vsbaseline.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
