#!/bin/bash

qsub -t 0-16   -l walltime=24:00:00,mem=16GB -o ./qsub_out/wrapper_invoperator.out -e ./qsub_out/wrapper_invoperator.err -v myvar='/home/rf02/rezvan/test1/step_by_step/Step4_MNE_Noise_Covariance_AND_Inverse_Operator.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
