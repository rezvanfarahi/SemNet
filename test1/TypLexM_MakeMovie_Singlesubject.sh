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

outpath='/imaging/rf02/TypLexMEG/singlesubject_maps/subtract'  # output directory for images
#inpath='/imaging/olaf/MEG/GoNoGo/STC/GM'
datapath='/imaging/rf02/TypLexMEG'
echo ${datapath}
#in_fif_fname='semloc_ssstf_fft49_raw.fif'
#eog_event_fname='semloc_ssstf_fft49_eog-eve.fif'
#ecg_event_fname='semloc_ssstf_fft49_ecg-eve.fif'
#out_fif_fname='semloc_ssstf_fft49_clean_ecg_eog_raw.fif'
#ecg_proj_fname='semloc_ssstf_fft49_ecg-proj.fif'
#eog_proj_fname='semloc_ssstf_fft49_eog-proj.fif'

# /imaging/olaf/MEG/GoNoGo/STC/GM/GM_lex_wdspds_n18.stc
########################

### All Conditions
conds=(\
'meg10_0378/101209/' \
'meg10_0390/101214/' \
'meg11_0026/110223/' \
'meg11_0050/110307/' \
'meg11_0052/110307/' \
'meg11_0069/110315/' \
'meg11_0086/110322/' \
'meg11_0091/110328/' \
'meg11_0096/110404/' \
'meg11_0101/110411/' \
'meg11_0102/110411/' \
'meg11_0112/110505/' \
'meg11_0104/110412/' \
'meg11_0118/110509/', \
'meg11_0131/110519/' \
'meg11_0144/110602/' \
'meg11_0147/110603/'
)

nconds=${#conds[*]}
lastcond=`expr $nconds - 1`

for cc in `seq 0 ${lastcond}`
do
  meg=${conds[cc]}
  aa='firstMorphed_SemLoc_icaclean_Subtract_Source_Evoked_meanresampled_0_550_100ms_50overlap'  
  echo ${meg}
  inpath=${datapath}/${meg}
  outdir=${outpath}/${meg}
  mkdir -p ${outdir}
  cd ${inpath}
  pwd
 
  mne_make_movie \
    --subject fsaverage \
    --smooth 5 \
    --stcin ${inpath}/${aa} \
    --scaleby 10000000000 \
    --tmin 0 \
    --tmax 450 \
    --tstep 50 \
    --pickrange \
    --alpha 1.0 \
    --width 600 \
    --height 400 \
    --noscalebar \
    --nocomments \
    --jpg ${outdir}/${aa} \
    --fthresh 0e10  --fmid 0.07 --fmax 0.15 
done    