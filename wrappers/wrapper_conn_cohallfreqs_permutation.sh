#!/bin/bash

qsub -t 0-14 -l walltime=24:00:00,mem=16GB -o ./qsub_out/wrapper_conn_cohallfreqs_permutatation.out -e ./qsub_out/wrapper_conn_cohallfreqs_permutatation.err -v myvar='/home/rf02/rezvan/semnet/avg_and_stat/Step3_4_MNE_Connectivity_cohallfreqs_evenodd_Permutation.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
