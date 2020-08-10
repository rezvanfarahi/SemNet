#!/bin/bash

qsub -t 9-16   -l walltime=24:00:00,mem=32GB -o ./qsub_out/wrapper_circmap.out -e ./qsub_out/wrapper_circmap.err -v myvar='/home/rf02/rezvan/test1/step_by_step/Step8_MNE_ROI_connectivity_wavelet.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
