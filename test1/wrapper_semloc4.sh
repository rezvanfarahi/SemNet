#!/bin/bash

qsub -t 0-16  -l walltime=24:00:00,mem=4GB -o ./qsub_out/wrapper_semloc4.out -e ./qsub_out/wrapper_semloc4.err -v myvar='/home/rf02/rezvan/test1/TyLexM_SourceBand_power4.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
