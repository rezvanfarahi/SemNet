#!/bin/bash

qsub -t 0-16   -l walltime=24:00:00,mem=32GB -o ./qsub_out/wrapper_wbconn_nonword.out -e ./qsub_out/wrapper_wbconn_nonword.err -v myvar='/home/rf02/rezvan/test1/step_by_step/TypLex/Step8_1_MNE_WB_connectivity_ctf_nonWord.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
