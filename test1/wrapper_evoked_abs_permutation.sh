#!/bin/bash

qsub -t 0-11   -l walltime=24:00:00,mem=16GB -o ./qsub_out/wrapper_evoked_abs_permutation.out -e ./qsub_out/wrapper_evoked_abs_permutation.err -v myvar='/home/rf02/rezvan/test1/TypLexM_Evoked_Cluster_ab_Permutation.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
