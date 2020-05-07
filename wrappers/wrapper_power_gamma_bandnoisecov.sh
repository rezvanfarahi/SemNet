#!/bin/bash

qsub -t 0-18  -l walltime=24:00:00,mem=32GB -o ./qsub_out/wrapper_power_gamma_bandnoisecov.out -e ./qsub_out/wrapper_power_gamma_bandnoisecov.err -v myvar='/home/rf02/rezvan/semnet/Step7_SourceBand_power_bandnoisecov_gamma.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
