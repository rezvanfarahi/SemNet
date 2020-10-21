#!/bin/bash

sbatch --job-name=Analysis1 --array 0-19 --time=24:00:00 --mem=32G -o /imaging/rf02/semnet_git_queue/sbatch_allinone_power_permutation.out -e /imaging/rf02/semnet_git_queue/sbatch_allinone_power_permutation.err  /imaging/rf02/semnet_git/SemNet/semnet4semloc/wrappers/wrapper_sbatch.sh '/imaging/rf02/semnet_git/SemNet/test1/step_by_step/averaging_and_stat/Step2_3_MNE_TF_Cluster_Permutation_allfreqs.py'
