#!/bin/bash

qsub -t 0-2   -l walltime=24:00:00,mem=16GB -o ./qsub_out/wrapper_coh150_permutation.out -e ./qsub_out/wrapper_coh150_permutation.err -v myvar='/home/rf02/rezvan/test1/replications/permutation/TypLexM_SpectralConnectivity_Permutation_coh_150.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
