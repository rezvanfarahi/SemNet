#!/bin/bash

qsub -t 0-50   -l walltime=24:00:00,mem=40GB -o ./qsub_out/wrapper_permutation_snr1.out -e ./qsub_out/wrapper_permutation_snr1.err -v myvar='/home/rf02/rezvan/test1/parcellation/extended_simulation_labels_ideal_snr1.py'  /home/rf02/rezvan/test1/wrapper_batch.sh