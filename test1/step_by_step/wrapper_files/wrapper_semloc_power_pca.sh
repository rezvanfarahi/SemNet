#!/bin/bash

qsub -t 0-12  -l walltime=24:00:00,mem=16GB -o ./qsub_out/wrapper_semloc_power_pca_permutation.out -e ./qsub_out/wrapper_semloc_power_pca_permutation.err -v myvar='/home/rf02/rezvan/test1/step_by_step/Step2_3_MNE_TF_pca_Cluster_Permutation.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
