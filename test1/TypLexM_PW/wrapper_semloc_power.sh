#!/bin/bash

qsub -t 0-15  -l walltime=24:00:00,mem=4GB -o ./qsub_out/wrapper_semloc_power.out -e ./qsub_out/wrapper_semloc_power.err -v myvar='/home/rf02/rezvan/test1/TypLexM_PW/TyLexM_SourceBand_power.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
