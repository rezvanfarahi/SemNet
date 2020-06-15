#!/bin/bash

sbatch --job-name=Analysis1 --array 0-19 -l --time=24:00:00 --mem=32G -o ./qsub_out/wrapper_apply_inverse_evoked.out  /imaging/rf02/semnet_git/SemNet/semnet4semloc/wrappers/wrapper_batch.sh '/imaging/rf02/semnet_git/SemNet/semnet4semloc/Step6_1_MNE_Applyinverse_evoked.py'
