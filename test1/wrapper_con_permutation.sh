#!/bin/bash

qsub -t 0-10   -l walltime=24:00:00,mem=16GB -o ./qsub_out/wrapper_con_permutation.out -e ./qsub_out/wrapper_con_permutation.err -v myvar='/home/rf02/rezvan/test1/TypLexM_SpectralConnectivity_Permutation1.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
