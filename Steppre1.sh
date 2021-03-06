#!/bin/bash



# Your variables
#before running script make sure you've typed:
# freesurfer_5.3.0
# mne_setup_2.7.3_64bit
# setenv SUBJECTS_DIR /imaging/rf02/Semnet/
datapath='/imaging/rf02/Semnet'       # root directory for your MEG data

# MEG IDs (your directory structure may differ)
subj_pre=(\
        'meg16_0030' \
        'meg16_0032' \
        'meg16_0033' \
	'meg16_0034' \
        'meg16_0035' \
        'meg16_0039' \
	'meg16_0041' \
        'meg16_0042' \
        'meg16_0045' \
	'meg16_0047' \
        'meg16_0052' \
	'meg16_0056' \
        'meg16_0069' \
        'meg16_0070' \
	'meg16_0071' \
        'meg16_0072' \
        'meg16_0073' \
	'meg16_0075' 
        )

# MEG subdirectories (your directory structure may differ)      
subj_dir=(\
         '160216' \
         '160218' \
         '160218' \
	 '160219' \
         '160222' \
	 '160225' \
	 '160226' \
	 '160229' \
	 '160303' \
	 '160304' \
	 '160310' \
	 '160314' \
	 '160405' \
	 '160407' \
	 '160407' \
	 '160408' \
	 '160411' \
	 '160411' 
        )

#condition names as used in file names
conds=('block_LD1_raw' 'block_LD2_raw' 'block_fruit_raw' 'block_milk_raw' 'block_odour_raw')
        

# Processing:

nsubjects=${#subj_pre[*]}
lastsubj=`expr $nsubjects - 1`

nconds=${#conds[*]}
lastcond=`expr $nconds - 1`


for m in `seq 0 ${lastsubj}`
do
  echo " "
  echo " Fixing electrodes for SUBJECT  ${subjects[m]}"
  echo " "
        for c in `seq 0 ${lastcond}`
        do
        
                mne_check_eeg_locations \
                        --file ${datapath}/${subj_pre[m]}/${subj_dir[m]}/${conds[c]}.fif \
                        --fix                           
                
        done # conditions

done # subjects
