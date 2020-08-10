#!/bin/bash

qsub -t 0-8   -l walltime=24:00:00,mem=32GB -o ./qsub_out/wrapper_circmap.out -e ./qsub_out/wrapper_circmap.err -v myvar='/home/rf02/rezvan/test1/TypLexM_circlemap_connectivity_wavelet2.py'  /home/rf02/rezvan/test1/wrapper_batch.sh
