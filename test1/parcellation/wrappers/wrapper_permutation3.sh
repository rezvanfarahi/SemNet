#!/bin/bash

qsub -t 10-36   -l walltime=24:00:00,mem=30GB -o ./qsub_out/wrapper_permutation3.out -e ./qsub_out/wrapper_permutation3.err -v myvar='/home/rf02/rezvan/test1/parcellation/extended_simulation_labels_ideal_fast_3seeds.py'  /home/rf02/rezvan/test1/wrapper_batch.sh