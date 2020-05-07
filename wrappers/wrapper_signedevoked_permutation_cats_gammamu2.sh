#!/bin/bash

qsub -t 0-5 -l walltime=24:00:00,mem=32GB -o ./qsub_out/wrapper_signedevoked_permutation_cats_gammamu2.out -e ./qsub_out/wrapper_signedevoked_permutation_cats_gammamu2.err -v myvar='/home/rf02/rezvan/semnet/avg_and_stat/Step3_1_MNE_signedEvoked_and_Cluster_Permutation_cats_gammamu2.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
