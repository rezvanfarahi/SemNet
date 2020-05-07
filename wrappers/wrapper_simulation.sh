#!/bin/bash

qsub -t 0-35   -l walltime=90:00:00,mem=40GB -o ./qsub_out/wrapper_simulation.out -e ./qsub_out/wrapper_simulation.err -v myvar='/home/rf02/rezvan/semnet/simulation_multivariate_coherence.py'  /home/rf02/rezvan/test1/wrapper_batch.sh