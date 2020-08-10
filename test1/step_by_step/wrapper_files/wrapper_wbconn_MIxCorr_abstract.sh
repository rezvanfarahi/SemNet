#!/bin/bash

qsub -t 0-16   -l walltime=24:00:00,mem=32GB -o ./qsub_out/wrapper_wbconn_MIxCorr_abstract.out -e ./qsub_out/wrapper_wbconn_MIxCorr_abstract.err -v myvar='/home/rf02/rezvan/test1/step_by_step/Step8_0_WB_Abstract_ctf_MI_connectivity.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
