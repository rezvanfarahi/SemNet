#!/bin/bash

sbatch --job-name=Analysis1 --array 0-19 --time=24:00:00 --mem=32G -o /imaging/rf02/semnet_git_queue/wrapper_apply_inverse_evoked.out -e /imaging/rf02/semnet_git_queue/wrapper_apply_inverse_evoked.err  /imaging/rf02/semnet_git/SemNet/semnet4semloc/wrappers/wrapper_sbatch.sh '/imaging/rf02/semnet_git/SemNet/semnet4semloc/Step6_1_MNE_Applyinverse_evoked.py'
