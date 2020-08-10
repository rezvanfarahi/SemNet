#!/bin/bash

qsub -t 0-16 -l walltime=24:00:00,mem=16GB -o ./qsub_out/wrapper_semloc_applyinverse1.out -e ./qsub_out/wrapper_semloc_applyinverse1.err -v myvar='/home/rf02/rezvan/test1/TypLexM_Applyinverse_epochs.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
