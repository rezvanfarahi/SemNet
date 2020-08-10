#!/bin/bash

qsub -t 0-12   -l walltime=24:00:00,mem=16GB -o ./qsub_out/wrapper_powtheta_permutation.out -e ./qsub_out/wrapper_powtheta_permutation.err -v myvar='/home/rf02/rezvan/test1/replications/permutation/TypLexM_SourceBand_Power_Permutation_theta.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
