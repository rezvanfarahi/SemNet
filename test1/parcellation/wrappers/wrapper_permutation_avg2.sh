#!/bin/bash

qsub -t 36-39   -l walltime=48:00:00,mem=40GB -o ./qsub_out/wrapper_permutation_avg2.out -e ./qsub_out/wrapper_permutation_avg2.err -v myvar='/home/rf02/rezvan/test1/parcellation/extended_simulation_labels_ideal_fast_avgspc2.py'  /home/rf02/rezvan/test1/wrapper_batch.sh