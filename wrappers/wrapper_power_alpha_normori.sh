#!/bin/bash

qsub -t 0-18  -l walltime=24:00:00,mem=32GB -o ./qsub_out/wrapper_power_alpha_normori.out -e ./qsub_out/wrapper_power_alpha_normori.err -v myvar='/home/rf02/rezvan/semnet/Step7_SourceBand_power_normori_alpha.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
