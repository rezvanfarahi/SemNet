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

outpath='/imaging/rf02/Semnet' #'/imaging/rf02/TypLexMEG/icaanalysis_results/jpg/permutation/connectivity/WB_hubs/to_use' #connectivity/WB_spokes/spectral'  # output directory for images
#inpath='/imaging/olaf/MEG/GoNoGo/STC/GM'
inpath='/mridata/cbu' #connectivity/WB_spokes/spectral/older'


# /imaging/olaf/MEG/GoNoGo/STC/GM/GM_lex_wdspds_n18.stc
########################
#'CBU160291_MR16002/20160415_085130/Series_002_xpace_MPRAGE_32chn/1.3.12.2.1107.5.2.43.67035.2016041509050189032922054.dcm' \
#'CBU140514_MR13008/20140609_154527/Series_003_MP_Rage_1x1x1/1.3.12.2.1107.5.2.32.35119.2014060915564950744748817.dcm' \
#'CBU160160_MR16001E/20160307_182629/Series_003_Accerelarated_Sctructural_Scan/1.3.12.2.1107.5.2.43.67035.2016030718440524425349752.dcm' \
#'CBU111174_MR11021/20111121_170544/Series_003_CBU_MPRAGE_32chn/1.3.12.2.1107.5.2.32.35119.2011112117124442561106614.dcm' \
#'CBU120947_MR12016/20120911_165153/Series_003_CBU_MPRAGE_32chn/1.3.12.2.1107.5.2.32.35119.2012091117121488427278414.dcm' \
#'CBU160144_MEG_STRUCTURALS/20160303_090816/Series_005_CBU_MPRAGE_32chn/1.3.12.2.1107.5.2.43.67035.2016030309143751427822748.dcm' \
#'CBU160310_MR16006/20160420_135651/Series_005_CBU_MPRAGE_32chn/1.3.12.2.1107.5.2.43.67035.2016042014164096737098.dcm' \
#'CBU131134_MR13018/20131217_113153/Series_003_CBU_MPRAGE_32chn/1.3.12.2.1107.5.2.32.35119.201312171201199932399939.dcm' \
#'CBU150562_MR13015/20151027_131116/Series_005_CBU_MPRAGE_32chn/1.3.12.2.1107.5.2.43.67035.2015102713183798492152591.dcm' \
#'CBU150064_METHODS/20150313_130602/Series_013_CBU_MPRAGE_32chn/1.3.12.2.1107.5.2.43.67035.201503131317335678306191.dcm' \
#'CBU160390_MEG_STRUCTURALS/20160510_111509/Series_005_CBU_MPRAGE_32chn/1.3.12.2.1107.5.2.43.67035.201605101122128660507637.dcm' \
#'CBU160208_MR13015/20160321_183603/Series_005_CBU_MPRAGE_32chn/1.3.12.2.1107.5.2.43.67035.201603211842289001859931.dcm' \
#'CBU111244_MR11021/20111213_170844/Series_008_CBU_MPRAGE_32chn/1.3.12.2.1107.5.2.32.35119.201112131834289163361206.dcm' \
#'CBU140887_MR12019/20140925_150421/Series_003_CBU_MPRAGE_32chn/1.3.12.2.1107.5.2.32.35119.201409251520407320027864.dcm' \
#'CBU160235_MR15007/20160331_164228/Series_005_CBU_MPRAGE_32chn/1.3.12.2.1107.5.2.43.67035.201603311710398371695109.dcm' \
#'CBU160385_CBU160385/20160509_120122/Series_005_CBU_MPRAGE_32chn/1.3.12.2.1107.5.2.43.67035.201605091207259343970744.dcm' \
#'CBU150037_MR15002/20150227_160047/Series_005_CBU_MPRAGE_64chn/1.3.12.2.1107.5.2.43.67035.2015022716090196934323359.dcm' \
#'CBU160304_MEG_STRUCTURAL/20160419_153358/Series_005_CBU_MPRAGE_32chn/1.3.12.2.1107.5.2.43.67035.201604191540237534714497.dcm' \
#'CBU160309_MEGSTRUCTURALS/20160420_132907/Series_005_CBU_MPRAGE_32chn/1.3.12.2.1107.5.2.43.67035.201604201335299124733994.dcm' \
#'CBU160373_MR16007/20160506_100711/Series_005_CBU_MPRAGE_32chn/1.3.12.2.1107.5.2.43.67035.201605061026242311106003.dcm' \
#'CBU150271_MR15008/20150605_163823/Series_005_CBU_MPRAGE_32chn/1.3.12.2.1107.5.2.43.67035.201506051648273735423166.dcm' \

### All Conditions
in_subs=(\

'CBU160610_MEG_STRUCTURALS/20160708_170933/Series_005_CBU_MPRAGE_32chn/1.3.12.2.1107.5.2.43.67035.2016070817165505401365.dcm' \
'CBU160646_MEG_STRUCURALS/20160721_172406/Series_005_CBU_MPRAGE_32chn/1.3.12.2.1107.5.2.43.67035.20160721173239880996558.dcm' \
'CBU160604_MR16007/20160707_132906/Series_005_CBU_MPRAGE_32chn/1.3.12.2.1107.5.2.43.67035.20160707135220781146242.dcm' \

)

out_subs=(\
#	    '/MRI_meg16_0030/mri/orig/001.mgz' \
#            '/MRI_meg16_0032/mri/orig/001.mgz' \
#            '/MRI_meg16_0033/mri/orig/001.mgz' \
#            '/MRI_meg16_0034/mri/orig/001.mgz' \
#            '/MRI_meg16_0035/mri/orig/001.mgz' \
#            '/MRI_meg16_0039/mri/orig/001.mgz' \
#            '/MRI_meg16_0041/mri/orig/001.mgz' \
#            '/MRI_meg16_0042/mri/orig/001.mgz' \
#            '/MRI_meg16_0045/mri/orig/001.mgz' \
#            '/MRI_meg16_0047/mri/orig/001.mgz' \
#            '/MRI_meg16_0052/mri/orig/001.mgz' \
#            '/MRI_meg16_0056/mri/orig/001.mgz' \
#            '/MRI_meg16_0069/mri/orig/001.mgz' \ 
#            '/MRI_meg16_0070/mri/orig/001.mgz' \
#            '/MRI_meg16_0072/mri/orig/001.mgz' \
#            '/MRI_meg16_0073/mri/orig/001.mgz' \
#            '/MRI_meg16_0075/mri/orig/001.mgz' \
#            '/MRI_meg16_0078/mri/orig/001.mgz' \
#            '/MRI_meg16_0082/mri/orig/001.mgz' \
#            '/MRI_meg16_0086/mri/orig/001.mgz' \
#            '/MRI_meg16_0097/mri/orig/001.mgz' \
	    '/MRI_meg16_0122/mri/orig/001.mgz' \
            '/MRI_meg16_0123/mri/orig/001.mgz' \
            '/MRI_meg16_0125/mri/orig/001.mgz' \
)

nconds=${#out_subs[*]}
lastcond=`expr $nconds - 1`

for cc in `seq 0 ${lastcond}`
do
  echo " Condition  ${cc}" 

  mri_convert ${inpath}/${in_subs[cc]} ${outpath}/${out_subs[cc]}
  
done    

