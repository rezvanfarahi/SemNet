#!/bin/bash

qsub -t 0-15  -l walltime=24:00:00,mem=4GB -o ./qsub_out/wrapper_semloc_power2.out -e ./qsub_out/wrapper_semloc_power2.err -v myvar='/home/rf02/rezvan/test1/TypLexM_PW/TyLexM_SourceBand_power2.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
