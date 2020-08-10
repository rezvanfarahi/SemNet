#!/bin/bash

qsub -t 0-16  -l walltime=24:00:00,mem=16GB -o ./qsub_out/wrapper_power_gamma.out -e ./qsub_out/wrapper_power_gamma.err -v myvar='/home/rf02/rezvan/test1/step_by_step/TypLex/Step7_4_MNE_SourceBand_power_gamma.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
