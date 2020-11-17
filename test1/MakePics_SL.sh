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

outpath='/imaging/rf02/TypLexMEG/icaanalysis_results/jpg/permutation/power/masked_ROIs_oct2020/20ms_wins/' #'/imaging/rf02/TypLexMEG/icaanalysis_results/jpg/permutation/connectivity/WB_hubs/to_use' #connectivity/WB_spokes/spectral'  # output directory for images
#inpath='/imaging/olaf/MEG/GoNoGo/STC/GM'
inpath='/imaging/rf02/TypLexMEG/icaanalysis_results/stc/permutation/power/masked_ROIs_oct2020/20ms_wins/' #'/imaging/rf02/TypLexMEG/icaanalysis_results/stc/UVttest/connectivity/WB_hubs/newp' #connectivity/WB_spokes/spectral/older'


# /imaging/olaf/MEG/GoNoGo/STC/GM/GM_lex_wdspds_n18.stc
########################

### All Conditions
conds=(\

##typlex power
'Per_sw_clusp05_p0424_SL_theta_ica_Subtract_Power2_ratio_normori_eq_50_550_20ms_sx_ms1' \

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
    --tmax 2450 \
    --tstep 50 \
    --pickrange \
    --alpha 1.0 \
    --width 600 \
    --height 400 \
    --noscalebar \
    --nocomments \
    --view  lat \
    --tif ${outpath}/${conds[cc]} \
    --fthresh 0e10  --fmid 0.01e10 --fmax 0.02e10
done    

#--fthresh 0e10  --fmid 0.01e10 --fmax 0.02e10
#--fthresh 0e10  --fmid 2e10 --fmax 4e10