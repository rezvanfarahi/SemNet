#!/bin/bash

qsub -t 0-2   -l walltime=24:00:00,mem=16GB -o ./qsub_out/wrapper_ppc250_permutation.out -e ./qsub_out/wrapper_ppc250_permutation.err -v myvar='/home/rf02/rezvan/test1/replications/permutation/TypLexM_SpectralConnectivity_Permutation_PPC_250.py'  /home/rf02/rezvan/test1/wrapper_batch.sh