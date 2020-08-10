#!/bin/bash

qsub -t 0-16   -l walltime=24:00:00,mem=32GB -o ./qsub_out/wrapper_wbconn_word.out -e ./qsub_out/wrapper_wbconn_word.err -v myvar='/home/rf02/rezvan/test1/step_by_step/TypLex/Step8_2_MNE_WB_connectivity_ctf_Word.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
