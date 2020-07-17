#!/bin/bash

sbatch --job-name=Permutation2 --array 3-4 --time=24:00:00 --mincpus=4 --mem-per-cpu=8G -o /imaging/rf02/semnet_git_queue/wrapper_rmanova_permutation.out -e /imaging/rf02/semnet_git_queue/wrapper_rmanova_permutation.err /imaging/rf02/semnet_git/SemNet/semnet4semloc/wrappers/wrapper_sbatch.sh '/imaging/rf02/semnet_git/SemNet/semnet4semloc/avg_and_stat/Step3_2_MNE_Evoked_Cluster_Permutation_rmanova.py' 
