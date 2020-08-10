#!/bin/bash

qsub -t 0-35   -l walltime=48:00:00,mem=40GB -o ./qsub_out/wrapper_permutation_avg3.out -e ./qsub_out/wrapper_permutation_avg3.err -v myvar='/home/rf02/rezvan/test1/parcellation/extended_simulation_labels_ideal_fast_avgspc3.py'  /home/rf02/rezvan/test1/wrapper_batch.sh