#!/bin/bash

qsub -t 0-24 -l walltime=24:00:00,mem=16GB -o ./qsub_out/wrapper_ica.out -e ./qsub_out/wrapper_ica.err -v myvar='/home/rf02/rezvan/semnet/Step3_MNE_ICA_par.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
