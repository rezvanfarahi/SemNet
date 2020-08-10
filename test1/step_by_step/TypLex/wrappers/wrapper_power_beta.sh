#!/bin/bash

qsub -t 0-16  -l walltime=24:00:00,mem=16GB -o ./qsub_out/wrapper_power_beta.out -e ./qsub_out/wrapper_power_beta.err -v myvar='/home/rf02/rezvan/test1/step_by_step/TypLex/Step7_3_MNE_SourceBand_power_beta.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
