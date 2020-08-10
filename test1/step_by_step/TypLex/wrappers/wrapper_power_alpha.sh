#!/bin/bash

qsub -t 0-16  -l walltime=24:00:00,mem=16GB -o ./qsub_out/wrapper_power_alpha.out -e ./qsub_out/wrapper_power_alpha.err -v myvar='/home/rf02/rezvan/test1/step_by_step/TypLex/Step7_2_MNE_SourceBand_power_alpha.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
