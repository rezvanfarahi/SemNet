#!/bin/bash

qsub -t 30-46   -l walltime=24:00:00,mem=30GB -o ./qsub_out/wrapper_permutation_extra.out -e ./qsub_out/wrapper_permutation_extra.err -v myvar='/home/rf02/rezvan/test1/parcellation/extended_simulation_labels_ideal_fast_15seeds.py'  /home/rf02/rezvan/test1/wrapper_batch.sh