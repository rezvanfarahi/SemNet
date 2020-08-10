#!/bin/bash

qsub -t 0-16 -l walltime=24:00:00,mem=16GB -o ./qsub_out/wrapper_applyinverse_epochs.out -e ./qsub_out/wrapper_applyinverse_epochs.err -v myvar='/home/rf02/rezvan/test1/step_by_step/TypLex/Step6_MNE_Applyinverse_epochs.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
