#!/bin/sh
#

## Your variables:

# datapath='/imaging/ec02/MEG/data/'    # root directory for your MEG data
datapath='/imaging/rf02/Semnet'

MRIpath='/imaging/rf02/Semnet'    # where your MRI subdirectories are

# subjects names used for MRI data
subjects=(\
            'MRI_meg16_0030' \
            'MRI_meg16_0032' \
            'MRI_meg16_0033' \
            'MRI_meg16_0034' \
            'MRI_meg16_0035' \
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

# MEG IDs (your directory structure may differ)
subj_pre=(\
	    'meg16_0030' \
            'meg16_0032' \
            'meg16_0033' \
            'meg16_0034' \
            'meg16_0035' \
            'meg16_0042' \
            'meg16_0045' \
            'meg16_0047' \
            'meg16_0052' \
            'meg16_0056' \
            'meg16_0069' \
 
            'meg16_0070' \
            'meg16_0072' \
            'meg16_0073' \
            'meg16_0075' \
            'meg16_0078' \
            'meg16_0082' \
            'meg16_0086' \
            'meg16_0097' \
            'meg16_0122' \
            'meg16_0123' \
            'meg16_0125' \
)
	 

# MEG subdirectories (your directory structure may differ)      
subj_dir=(\
	    '160216' \
            '160218' \
            '160218' \
            '160219' \
            '160222' \
            '160229' \
            '160303' \
            '160304' \
            '160310' \
            '160314' \
            '160405' \
            '160407' \
            '160408' \
            '160411' \
            '160411' \
            '160414' \
            '160418' \
            '160422' \
            '160512' \
            '160707' \
            '160708' \
            '160712' \
        )
        

## Processing:

nsubjects=${#subjects[*]}
lastsubj=`expr $nsubjects - 1`

for m in `seq 0 ${lastsubj}`
do
  echo " "
  echo " Computing forward solution for SUBJECT  ${subjects[m]}"
  echo " "

                           
  ## setup model 3 layers (EEG+MEG)
  
  mne_setup_forward_model --overwrite  --subject ${subjects[m]} --surf --ico 4
    
  mne_do_forward_solution \
                        --overwrite \
                        --subject ${subjects[m]} \
                        --mindist 5 \
                        --spacing ico-5 \
                        --bem ${MRIpath}/${subjects[m]}/bem/${subjects[m]}-5120-5120-5120-bem-sol.fif \
                        --src ${MRIpath}/${subjects[m]}/bem/${subjects[m]}-ico-5-src.fif \
                        --meas ${datapath}/${subj_pre[m]}/${subj_dir[m]}/block_LD1_raw.fif  \
                        --fwd ${datapath}/${subj_pre[m]}/${subj_dir[m]}/ico5_forward_5-3L-EMEG-fwd.fif




done # subjects


#--trans ${MRIpath}/${subjects[m]}/bem/${subjects[m]}-trans.fif