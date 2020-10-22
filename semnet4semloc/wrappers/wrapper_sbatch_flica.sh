#!/bin/bash

myvar=$1
echo "SLURM_ARRAY_TASK_ID " $SLURM_ARRAY_TASK_ID
echo $myvar

# myvar: python script to execute
# PBS_ARRAYID: subject index to be processed

#/imaging/local/software/anaconda/2.4.1/2/bin/python2.7 ${myvar} $SLURM_ARRAY_TASK_ID
/home/rf02/.conda/envs/profumo_flica/bin/python2.7 ${myvar} $SLURM_ARRAY_TASK_ID
