#!/bin/bash

qsub -t 0-18  -l walltime=24:00:00,mem=32GB -o ./qsub_out/wrapper_power_theta_normori.out -e ./qsub_out/wrapper_power_theta_normori.err -v myvar='/home/rf02/rezvan/semnet/Step7_SourceBand_power_normori_theta.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
