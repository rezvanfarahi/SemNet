#!/bin/bash

sbatch --job-name=Permutation1 --array 0-6 --time=24:00:00 --mem=32G -o /imaging/rf02/semnet_git_queue/wrapper_evoked_permutation.out -e /imaging/rf02/semnet_git_queue/wrapper_evoked_permutation.err /imaging/rf02/semnet_git/SemNet/semnet4semloc/wrappers/wrapper_sbatch.sh '/imaging/rf02/semnet_git/SemNet/semnet4semloc/avg_and_stat/Step3_1_MNE_Evoked_Cluster_Permutation.py' 
