#!/bin/bash

qsub -t 0-16  -l walltime=24:00:00,mem=16GB -o ./qsub_out/wrapper_semloc06.out -e ./qsub_out/wrapper_semloc06.err -v myvar='/home/rf02/rezvan/test1/replications/epochs/TypLexM_AbstractBaseline_SpectralConn_SourceSpace_MNE.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
