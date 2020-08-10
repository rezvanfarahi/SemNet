#!/bin/bash

qsub -t 0-16   -l walltime=24:00:00,mem=32GB -o ./qsub_out/wrapper_wbconn_temp_concrete_wband.out -e ./qsub_out/wrapper_wbconn_temp_concrete_wband.err -v myvar='/home/rf02/rezvan/test1/step_by_step/Step10_2_WB_Concrete_ctf_MI_connectivity_wband.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
