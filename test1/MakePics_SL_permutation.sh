#!/bin/sh
#

#before running script make sure you've typed:
# freesurfer_4.3.0
# mne_setup_2.7.3_64bit
# setenv SUBJECTS_DIR /imaging/rf02/TypLexMEG/


################
# PARAMETERS:

# Input/Output file path:

#path='/megdata/cbu/PATH/TO/YOUR/DATA'
#path='/group/erp/data/olaf.hauk/Others/Miozzo/data'
#MRIpath='/group/erp/data/olaf.hauk/Others/Miozzo/MRIs'

outpath='/imaging/rf02/TypLexMEG/output_stc_imageseries/SemLoc/permutationmaps'  # output directory for images
#inpath='/imaging/olaf/MEG/GoNoGo/STC/GM'
inpath='/imaging/rf02/TypLexMEG'

# /imaging/olaf/MEG/GoNoGo/STC/GM/GM_lex_wdspds_n18.stc
########################

### All Conditions
conds=(\
#'Permutation_ttest_pmap_SemLoc_Concrete_SourcePower_Logratio_m500_700_theta' \
#'Permutation_ttest_pmap_SemLoc_Abstract_SourcePower_Logratio_m500_700_theta' \
#'Permutation_ttest_pmap_SemLoc_Concrete_SourcePower_Logratio_m500_700_alpha' \
#'Permutation_ttest_pmap_SemLoc_Abstract_SourcePower_Logratio_m500_700_alpha' \
#'Permutation_ttest_pmap_SemLoc_Concrete_SourcePower_Logratio_m500_700_beta' \
#'Permutation_ttest_pmap_SemLoc_Abstract_SourcePower_Logratio_m500_700_beta' \
#'Permutation_ttest_pmap_SemLoc_Concrete_SourcePower_Logratio_m500_700_gamma' \
#'Permutation_ttest_pmap_SemLoc_Abstract_SourcePower_Logratio_m500_700_gamma' \
'Permutation_clusterp005_p05_15subj_SemLoc_icaclean_ConcreteAbstract_SourcePower_ratio_0_600_theta'\
)

nconds=${#conds[*]}
lastcond=`expr $nconds - 1`
#pmap medial
for cc in `seq 0 ${lastcond}`
do
  echo " Condition  ${cc}" 

  mne_make_movie \
    --subject avgsubject \
    --smooth 5 \
    --stcin ${inpath}/${conds[cc]} \
    --tmin 0 \
    --tmax 1 \
    --tstep 1 \
    --pickrange \
    --alpha 1.0 \
    --width 600 \
    --height 400 \
    --nocomments \
    --noscalebar \
    --jpg ${outpath}/${conds[cc]} \
    --fthresh 0e10  --fmid 100e10  --fmax 120e10 
done    

#pmap lateral
for cc in `seq 0 ${lastcond}`
do
  echo " Condition  ${cc}" 

  mne_make_movie \
    --subject avgsubject \
    --smooth 5 \
    --stcin ${inpath}/${conds[cc]} \
    --tmin 0 \
    --tmax 6 \
    --tstep 1 \
    --pickrange \
    --alpha 1.0 \
    --width 600 \
    --height 400 \
    --nocomments \
    --noscalebar \
    --view med \
    --jpg ${outpath}/${conds[cc]} \
    --fthresh 0.95e10  --fmid 0.975e10  --fmax 1e10 
done    

conds=(\
'Permutation_ttest_tmap_SemLoc_Concrete_SourcePower_Logratio_m500_700_theta' \
'Permutation_ttest_tmap_SemLoc_Abstract_SourcePower_Logratio_m500_700_theta' \
'Permutation_ttest_tmap_SemLoc_Concrete_SourcePower_Logratio_m500_700_alpha' \
'Permutation_ttest_tmap_SemLoc_Abstract_SourcePower_Logratio_m500_700_alpha' \
'Permutation_ttest_tmap_SemLoc_Concrete_SourcePower_Logratio_m500_700_beta' \
'Permutation_ttest_tmap_SemLoc_Abstract_SourcePower_Logratio_m500_700_beta' \
'Permutation_ttest_tmap_SemLoc_Concrete_SourcePower_Logratio_m500_700_gamma' \
'Permutation_ttest_tmap_SemLoc_Abstract_SourcePower_Logratio_m500_700_gamma' \
)

nconds=${#conds[*]}
lastcond=`expr $nconds - 1`

#tmap medial
for cc in `seq 0 ${lastcond}`
do
  echo " Condition  ${cc}" 

  mne_make_movie \
    --subject avgsubject \
    --smooth 5 \
    --stcin ${inpath}/${conds[cc]} \
    --tmin 0 \
    --tmax 6 \
    --tstep 1 \
    --pickrange \
    --alpha 1.0 \
    --width 600 \
    --height 400 \
    --nocomments \
    --noscalebar \
    --view med \
    --jpg ${outpath}/${conds[cc]} \
    --fthresh 5e10  --fmid 7e10  --fmax 9e10 
done    

#tmap lateral
for cc in `seq 0 ${lastcond}`
do
  echo " Condition  ${cc}" 

  mne_make_movie \
    --subject avgsubject \
    --smooth 5 \
    --stcin ${inpath}/${conds[cc]} \
    --tmin 0 \
    --tmax 6 \
    --tstep 1 \
    --pickrange \
    --alpha 1.0 \
    --width 600 \
    --height 400 \
    --nocomments \
    --noscalebar \
    --view lat \
    --jpg ${outpath}/${conds[cc]} \
    --fthresh 5e10  --fmid 7e10  --fmax 9e10 
done    

