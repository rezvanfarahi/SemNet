#!/bin/bash

qsub -t 17 -l walltime=24:00:00,mem=32GB -o ./qsub_out/wrapper_apply_inverse_epochs_spec_conn2.out -e ./qsub_out/wrapper_apply_inverse_epochs_spec_conn2.err -v myvar='/home/rf02/rezvan/semnet/Step8_0_WB_applyinv_epochs_Coh_meanctf_connectivity.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
