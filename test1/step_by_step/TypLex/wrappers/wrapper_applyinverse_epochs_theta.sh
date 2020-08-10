#!/bin/bash

qsub -t 0-16 -l walltime=24:00:00,mem=32GB -o ./qsub_out/wrapper_applyinverse_epochs_theta.out -e ./qsub_out/wrapper_applyinverse_epochs_theta.err -v myvar='/home/rf02/rezvan/test1/step_by_step/TypLex/Step9_1_MNE_Applyinverse_epochs_theta.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
