#!/bin/sh
#
#before running script make sure you've typed:
# freesurfer_5.3.0
# mne_setup_2.7.3_64bit
# setenv SUBJECTS_DIR /imaging/rf02/TypLexMEG/

## Your variables:

# datapath='/imaging/ec02/MEG/data/'    # root directory for your MEG data
datapath='/home/rf02/rezvan/TypLexMEG'

MRIpath='/home/rf02/rezvan/TypLexMEG'    # where your MRI subdirectories are

# subjects names used for MRI data
subjects=(\
         'SS1_Meg10_0378' \
	'SS2_Meg10_0390' \
	'SS3_Meg11_0026' \
	'SS4_Meg11_0050' \
	'SS5_Meg11_0052' \
	'SS6_Meg11_0069' \
	'SS9_Meg11_0086' \
	'SS10_Meg11_0091' \
	'SS12_Meg11_0096' \
	'SS13_Meg11_0101' \
	'SS14_Meg11_0102' \
	'SS15_Meg11_0112' \
	'SS16_Meg11_0104' \
	'SS18_Meg11_0118' \
	'SS19_Meg11_0131' \
	'SS20_Meg11_0144' \
	'SS21_Meg11_0147' \
        )


## Processing:

make_average_subject \
		   --out avgsubj \
		   --subjects ${subjects[0]} ${subjects[1]} ${subjects[2]} ${subjects[3]} ${subjects[4]} \ ${subjects[5]} ${subjects[6]} ${subjects[7]} ${subjects[8]} ${subjects[9]} ${subjects[10]} \
${subjects[11]} ${subjects[12]} ${subjects[13]} ${subjects[14]} ${subjects[15]} ${subjects[16]} 

tkregister2 --fstal --s avgsubj --mgz 

