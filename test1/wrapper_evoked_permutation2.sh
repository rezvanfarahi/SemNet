#!/bin/bash

qsub -t 0-16   -l walltime=24:00:00,mem=16GB -o ./qsub_out/wrapper_evoked_permutation2.out -e ./qsub_out/wrapper_evoked_permutation2.err -v myvar='/home/rf02/rezvan/test1/TypLexM_Evoked_Cluster_Permutation2.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
