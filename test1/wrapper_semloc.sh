#!/bin/bash

qsub -t 0-1   -l walltime=24:00:00,mem=4GB -o ./qsub_out/wrapper_semloc.out -e ./qsub_out/wrapper_semloc.err -v myvar='/home/rf02/rezvan/test1/TyLexM_SourceBand_power.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
