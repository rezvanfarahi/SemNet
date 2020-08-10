#!/bin/bash

qsub -t 0-40   -l walltime=24:00:00,mem=30GB -o ./qsub_out/wrapper_permutation_null.out -e ./qsub_out/wrapper_permutation_null.err -v myvar='/home/rf02/rezvan/test1/parcellation/extended_simulation_labels_ideal_null.py'  /home/rf02/rezvan/test1/wrapper_batch.sh