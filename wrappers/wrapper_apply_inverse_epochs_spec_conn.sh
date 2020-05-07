#!/bin/bash

qsub -t 0-18 -l walltime=124:00:00,mem=64GB -o ./qsub_out/wrapper_apply_inverse_epochs_spec_conn.out -e ./qsub_out/wrapper_apply_inverse_epochs_spec_conn.err -v myvar='/home/rf02/rezvan/semnet/Step8_0_WB_applyinv_epochs_Coh_meanctf_allverts_connectivity.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
