#!/bin/sh
#

#before running script make sure you've typed:
# freesurfer_5.3.0
# mne_setup_2.7.3_64bit
# setenv SUBJECTS_DIR /imaging/rf02/TypLexMEG/
#NOTE! use TypLexMEG folder not Semnet (for fsaverage)

################
# PARAMETERS:

# Input/Output file path:

#path='/megdata/cbu/PATH/TO/YOUR/DATA'
#path='/group/erp/data/olaf.hauk/Others/Miozzo/data'
#MRIpath='/group/erp/data/olaf.hauk/Others/Miozzo/MRIs'

outpath='/imaging/rf02/Semnet/semnet4semloc/jpg/permutation/evoked/glm/' 
#'/imaging/rf02/Semnet/semnet4semloc/jpg/uvttest/evoked/glm' #'/imaging/rf02/Semnet/jpg/uvttest/evoked/pnt1_30/newregsigned/cats/Hand_Hear' #'/imaging/rf02/Semnet/jpg/uvttest/evoked/bands/newregsigned/mugamma' #'/imaging/rf02/Semnet/jpg/permutation/evoked/bands/newregsigned/mugamma' #'/imaging/rf02/Semnet/jpg/uvttest/evoked/pnt1_30/newregsigned/cats' #'/imaging/rf02/Semnet/jpg/uvttest/evoked/pnt1_30/newregsigned/wpw' #'/imaging/rf02/Semnet/jpg/GrandAverage/evoked/pnt1_30/newregsigned/cats/Visual' #'/imaging/rf02/Semnet/jpg/uvttest/evoked/pnt1_30/newregsigned/cats'  # output directory for images
#inpath='/imaging/olaf/MEG/GoNoGo/STC/GM'
inpath='/imaging/rf02/Semnet/semnet4semloc/stc/permutation/evoked/glm/' 
#'/imaging/rf02/Semnet/semnet4semloc/stc/uvttest/evoked/glm' #'/imaging/rf02/Semnet/stc/uvttest/evoked/pnt1_30/wbwnewregsigned2/cats' #'/imaging/rf02/Semnet/stc/uvttest/evoked/bands/newregsigned/mugamma' #'/imaging/rf02/Semnet/stc/permutation/evoked/bands/newregsigned/mugamma' #'/imaging/rf02/Semnet/stc/uvttest/evoked/pnt1_30/newregsigned/cats' #'/imaging/rf02/Semnet/stc/uvttest/evoked/pnt1_30/newregsigned/wpw' #'/imaging/rf02/Semnet/stc/uvttest/evoked/pnt1_30/newregsigned/cats' #


# /imaging/olaf/MEG/GoNoGo/STC/GM/GM_lex_wdspds_n18.stc
########################

### All Conditions
conds=(\

'ClusPer_GLM_Evoked_sw_icomorphed_oldreg_clusterp05_p04519096180763847_53subj_pnt1_48ica_contrast_5' \
#'UVTtest_t_icomorphed_oldreg_53subj_GLM_50_450_100ms_1_48ica_contrast' \

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
    --scaleby 10000000000 \
    --tmin 50 \
    --tmax 400 \
    --tstep 100 \
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