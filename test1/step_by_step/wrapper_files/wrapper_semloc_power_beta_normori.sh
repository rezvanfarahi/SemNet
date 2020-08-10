#!/bin/bash

qsub -t 8-15  -l walltime=24:00:00,mem=32GB -o ./qsub_out/wrapper_semloc_power_beta_normori.out -e ./qsub_out/wrapper_semloc_power_beta_normori.err -v myvar='/home/rf02/rezvan/test1/step_by_step/Step6new_SourceBand_power_normori_beta.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
