#!/bin/bash

qsub -t 0-16 -l walltime=24:00:00,mem=32GB -o ./qsub_out/wrapper_applyinverse_epochsa.out -e ./qsub_out/wrapper_applyinverse_epochsa.err -v myvar='/home/rf02/rezvan/test1/step_by_step/Step10_1_MNE_Applyinverse_epochs_alpha.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
