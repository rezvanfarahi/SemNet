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

outpath='/imaging/rf02/Semnet/semnet4semloc/stc/uvttest/evoked/1_48/ConEmot/' #'/imaging/rf02/Semnet/jpg/uvttest/evoked/pnt1_30/newregsigned/cats/Hand_Hear' #'/imaging/rf02/Semnet/jpg/uvttest/evoked/bands/newregsigned/mugamma' #'/imaging/rf02/Semnet/jpg/permutation/evoked/bands/newregsigned/mugamma' #'/imaging/rf02/Semnet/jpg/uvttest/evoked/pnt1_30/newregsigned/cats' #'/imaging/rf02/Semnet/jpg/uvttest/evoked/pnt1_30/newregsigned/wpw' #'/imaging/rf02/Semnet/jpg/GrandAverage/evoked/pnt1_30/newregsigned/cats/Visual' #'/imaging/rf02/Semnet/jpg/uvttest/evoked/pnt1_30/newregsigned/cats'  # output directory for images
#inpath='/imaging/olaf/MEG/GoNoGo/STC/GM'
inpath='/imaging/rf02/Semnet/semnet4semloc/stc/uvttest/evoked/1_48' #'/imaging/rf02/Semnet/stc/uvttest/evoked/pnt1_30/wbwnewregsigned2/cats' #'/imaging/rf02/Semnet/stc/uvttest/evoked/bands/newregsigned/mugamma' #'/imaging/rf02/Semnet/stc/permutation/evoked/bands/newregsigned/mugamma' #'/imaging/rf02/Semnet/stc/uvttest/evoked/pnt1_30/newregsigned/cats' #'/imaging/rf02/Semnet/stc/uvttest/evoked/pnt1_30/newregsigned/wpw' #'/imaging/rf02/Semnet/stc/uvttest/evoked/pnt1_30/newregsigned/cats' #


# /imaging/olaf/MEG/GoNoGo/STC/GM/GM_lex_wdspds_n18.stc
########################

### All Conditions
conds=(\


'UVTtest_icomorphed_oldreg_19subj_SemDec_m300_600_100ms_SL_1_48ica_Concrete_Emotional' \
'UVTtest_icomorphed_oldreg_19subj_SemDec_m300_600_100ms_SL_1_48ica_Concrete_Emotional' \


##perm wpw
#'ClusPer_abs_sw_icomorphed_newreg_clusterp05_p05_19subj_SemDec_pnt1_30ica_Hand_Pwordc_maxstep5' \
#'ClusPer_abs_sw_icomorphed_newreg_clusterp05_p05_19subj_SemDec_pnt1_30ica_Visual_Pwordc_maxstep5' \

#'ClusPer_abs_sw_icomorphed_newreg_clusterp05_p05_19subj_SemDec_pnt1_30ica_Visual_Hand_maxstep5' \
#'ClusPer_sEvoked_sw_icomorphed_newreg_clusterp01_p05_19subj_SemDec_pnt1_30ica_Visual_Handgammamu_maxstep5' \


)

nconds=${#conds[*]}
lastcond=`expr $nconds - 1`

for cc in `seq 0 ${lastcond}`
do
  echo " Condition  ${cc}"
  echo " ${outpath}/${conds[cc]} " 

  mne_make_movie \
    --subject fsaverage \
    --smooth 5 \
    --stcin ${inpath}/${conds[cc]} \
    --tmin 50 \
    --tmax 500 \
    --tstep 100 \
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
