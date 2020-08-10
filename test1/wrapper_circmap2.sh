#!/bin/bash

qsub -t 9-16   -l walltime=24:00:00,mem=32GB -o ./qsub_out/wrapper_circmap2.out -e ./qsub_out/wrapper_circmap2.err -v myvar='/home/rf02/rezvan/test1/TypLexM_circlemap_connectivity_wavelet.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
