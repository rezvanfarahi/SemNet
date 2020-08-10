#!/bin/bash

qsub -t 0-15   -l walltime=24:00:00,mem=40GB -o ./qsub_out/wrapper_permutation15.out -e ./qsub_out/wrapper_permutation15.err -v myvar='/home/rf02/rezvan/test1/parcellation/extended_simulation_labels_ideal_fast_15seeds.py'  /home/rf02/rezvan/test1/wrapper_batch.sh