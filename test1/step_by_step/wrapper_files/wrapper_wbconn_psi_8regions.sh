#!/bin/bash

qsub -t 0-16   -l walltime=24:00:00,mem=32GB -o ./qsub_out/wrapper_wbconn_psi_8regions.out -e ./qsub_out/wrapper_wbconn_psi_8regions.err -v myvar='/home/rf02/rezvan/test1/step_by_step/Step8_MNE_WB2spokes_connectivity_PSI.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
