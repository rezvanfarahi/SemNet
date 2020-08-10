#!/bin/bash

qsub -t 0-16  -l walltime=24:00:00,mem=16GB -o ./qsub_out/wrapper_semloc01.out -e ./qsub_out/wrapper_semloc01.err -v myvar='/home/rf02/rezvan/test1/replications/epochs/TypLexM_Concrete_wband_Coherence_SourceSpace_MNE_150_350.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
