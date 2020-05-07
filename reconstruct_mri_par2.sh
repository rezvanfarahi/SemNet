#!/bin/bash
#PBS -q compute
#PBS -l walltime=24:00:00,mem=64GB
#PBS -o ~/rezvan/semnet/reconstruct_mri_bash/reconstruct_mri_autorecon2.out 
#PBS -e ~/rezvan/semnet/reconstruct_mri_bash/reconstruct_mri_autorecon2.err



source /imaging/local/software/freesurfer/5.3.0/x86_64/SetUpFreeSurfer.sh >> /dev/null
export SUBJECTS_DIR=/imaging/rf02/Semnet
echo $SUBJECTS_DIR

out_subs=(\
	    'MRI_meg16_0030' \
            'MRI_meg16_0032' \
            'MRI_meg16_0033' \
            'MRI_meg16_0034' \
            'MRI_meg16_0035' \
            'MRI_meg16_0039' \
            'MRI_meg16_0041' \
            'MRI_meg16_0042' \
            'MRI_meg16_0045' \
            'MRI_meg16_0047' \
            'MRI_meg16_0052' \
            'MRI_meg16_0056' \
	    'MRI_meg16_0069' \
 	    'MRI_meg16_0070' \
            'MRI_meg16_0072' \
            'MRI_meg16_0073' \
            'MRI_meg16_0075' \
            'MRI_meg16_0078' \
            'MRI_meg16_0082' \
            'MRI_meg16_0086' \
            'MRI_meg16_0097' \
            'MRI_meg16_0122' \
            'MRI_meg16_0123' \
            'MRI_meg16_0125' \
            
)

current_sub=${out_subs[$PBS_ARRAYID]}
echo "subject $current_sub" 
recon-all -subjid $current_sub -autorecon3
  
   

