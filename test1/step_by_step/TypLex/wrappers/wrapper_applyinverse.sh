#!/bin/bash

qsub -t 0-16   -l walltime=24:00:00,mem=16GB -o ./qsub_out/wrapper_applyinverse.out -e ./qsub_out/wrapper_applyinverse.err -v myvar='/home/rf02/rezvan/test1/step_by_step/TypLex/Step5_MNE_Applyinverse_evoked.py'  /home/rf02/rezvan/test1/wrapper_batch.sh