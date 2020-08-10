#!/bin/bash

qsub -t 0-16   -l walltime=24:00:00,mem=4GB -o ./qsub_out/wrapper_ecgeog_ssp.out -e ./qsub_out/wrapper_ecgeog_ssp.err -v myvar='/home/rf02/rezvan/test1/TypLexM_eogecg_ssp_parallel.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
