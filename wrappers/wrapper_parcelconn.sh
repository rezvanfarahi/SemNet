#!/bin/bash

qsub -t 0-19 -l walltime=24:00:00,mem=32GB -o ./qsub_out/wrapper_parcelconn.out -e ./qsub_out/wrapper_parcelconn.err -v myvar='/home/rf02/rezvan/semnet/Step9_parcellation_connectivity.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
