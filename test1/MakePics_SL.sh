#!/bin/sh
#

#before running script make sure you've typed:
# freesurfer_5.3.0
# mne_setup_2.7.3_64bit
# setenv SUBJECTS_DIR /imaging/rf02/TypLexMEG/


################
# PARAMETERS:

# Input/Output file path:

#path='/megdata/cbu/PATH/TO/YOUR/DATA'
#path='/group/erp/data/olaf.hauk/Others/Miozzo/data'
#MRIpath='/group/erp/data/olaf.hauk/Others/Miozzo/MRIs'

outpath='/imaging/rf02/TypLexMEG/icaanalysis_results/jpg/UVttest/power/20ms_wins' #'/imaging/rf02/TypLexMEG/icaanalysis_results/jpg/permutation/connectivity/WB_hubs/to_use' #connectivity/WB_spokes/spectral'  # output directory for images
#inpath='/imaging/olaf/MEG/GoNoGo/STC/GM'
inpath='/imaging/rf02/TypLexMEG/icaanalysis_results/stc/UVttest/power/20ms_wins' #'/imaging/rf02/TypLexMEG/icaanalysis_results/stc/UVttest/connectivity/WB_hubs/newp' #connectivity/WB_spokes/spectral/older'


# /imaging/olaf/MEG/GoNoGo/STC/GM/GM_lex_wdspds_n18.stc
########################

### All Conditions
conds=(\

##typlex power
'UVTtest_t_Power2_ratio_normori_eq_50_550_20ms_theta' \

)

nconds=${#conds[*]}
lastcond=`expr $nconds - 1`

for cc in `seq 0 ${lastcond}`
do
  echo " Condition  ${cc}" 

  mne_make_movie \
    --subject fsaverage \
    --smooth 5 \
    --stcin ${inpath}/${conds[cc]} \
    --scaleby 10000000000 \
    --tmin 50 \
    --tmax 550 \
    --tstep 20 \
    --pickrange \
    --alpha 1.0 \
    --width 600 \
    --height 400 \
    --noscalebar \
    --nocomments \
    --view  lat \
    --jpg ${outpath}/${conds[cc]} \
    --fthresh 0e10  --fmid 2e10 --fmax 4e10
done    

