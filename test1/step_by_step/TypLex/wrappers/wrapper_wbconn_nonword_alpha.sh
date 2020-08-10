#!/bin/bash

qsub -t 0-16   -l walltime=24:00:00,mem=32GB -o ./qsub_out/wrapper_wbconn_nonword_alpha.out -e ./qsub_out/wrapper_wbconn_nonword_alpha.err -v myvar='/home/rf02/rezvan/test1/step_by_step/TypLex/Step8_1_MNE_WB_MI_connectivity_ctf_nonWord_alpha.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
