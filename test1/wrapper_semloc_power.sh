#!/bin/bash

qsub -t 3  -l walltime=24:00:00,mem=2GB -o ./qsub_out/wrapper_semloc_power.out -e ./qsub_out/wrapper_semloc_power.err -v myvar='/home/rf02/rezvan/test1/TyLexM_SourceBand_power.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
 