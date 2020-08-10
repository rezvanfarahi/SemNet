#!/bin/bash

qsub -t 35-70   -l walltime=24:00:00,mem=40GB -o ./qsub_out/wrapper_permutation2.out -e ./qsub_out/wrapper_permutation2.err -v myvar='/home/rf02/rezvan/test1/parcellation/extended_simulation_labels_ideal_fast2.py'  /home/rf02/rezvan/test1/wrapper_batch.sh