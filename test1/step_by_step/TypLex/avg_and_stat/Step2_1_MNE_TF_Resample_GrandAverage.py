"""
=========================================================
TF representation of SemLoc for frequency bands in source labels
=========================================================

"""
# Authors: Alexandre Gramfort <gramfort@nmr.mgh.harvard.edu>
#
# License: BSD (3-clause)


print(__doc__)

import sys
sys.path.append('/imaging/local/software/python_packages/nibabel/1.3.0')
sys.path.append('/imaging/local/software/python_packages/pysurfer/v0.3.1')
sys.path.append('/imaging/local/software/python_packages/scikit-learn/v0.14.1')
sys.path.append('/imaging/local/software/mne_python/latest')

# for qsub
# (add ! here if needed) 
sys.path.append('/imaging/local/software/anaconda/latest/x86_64/bin/python')
sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
###

import numpy as np
import mne
import os


###############################################################################
out_path = '/imaging/rf02/TypLexMEG/icaanalysis_results/typlex/stc/GrandAverage/power/final/' # root directory for your MEG data
if not os.path.exists(out_path):
    os.makedirs(out_path)
subjects_dir = '/imaging/rf02/TypLexMEG/'    # where your MRI subdirectories are
# where event-files are
data_path = '/imaging/rf02/TypLexMEG/'    # where event files are

label_path = '/home/rf02/rezvan/TypLexMEG/createdlabels_SL/'

inv_path = '/imaging/rf02/TypLexMEG/'
inv_fname = 'InverseOperator_EMEG-inv.fif'

# get indices for subjects to be processed from command line input
# 
print sys.argv
subject_inds = []
for ss in sys.argv[1:]:
   subject_inds.append( int( ss ) )

subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
print "subject_inds:"
print subject_inds
print "No rejection"

list_all = ['meg10_0378/101209/', 
'meg10_0390/101214/',
'meg11_0026/110223/', 
'meg11_0050/110307/', 
'meg11_0052/110307/', 
'meg11_0069/110315/', 
'meg11_0086/110322/', 
'meg11_0091/110328/', 
'meg11_0096/110404/', 
'meg11_0101/110411/', 
'meg11_0102/110411/', 
'meg11_0112/110505/',
'meg11_0104/110412/',  
'meg11_0118/110509/', 
'meg11_0131/110519/', 
'meg11_0144/110602/', 
'meg11_0147/110603/', 


]

# subjects names used for MRI data
subjects=['SS1_Meg10_0378',
'SS2_Meg10_0390',
'SS3_Meg11_0026',
'SS4_Meg11_0050',
'SS5_Meg11_0052',
'SS6_Meg11_0069',
'SS9_Meg11_0086',
'SS10_Meg11_0091',
'SS12_Meg11_0096',
'SS13_Meg11_0101',
'SS14_Meg11_0102',
'SS15_Meg11_0112',
'SS16_Meg11_0104',
'SS18_Meg11_0118',
'SS19_Meg11_0131',
'SS20_Meg11_0144',
'SS21_Meg11_0147'
]

ll = []
for ss in subject_inds:
 ll.append(list_all[ss])

print "ll:"
print ll



# Labels/ROIs
#labellist = ['atlleft-lh']#, 'atlleftt-rh'
tmin, tmax = -0.5, 0.7
# frequency bands with variable number of cycles for wavelets
stim_delay = 0.034 # delay in s


tmin1=50
tstep1=100
stc_allc=range(len(subject_inds))
vertices_to = [np.arange(10242), np.arange(10242)]
stc_alla=range(len(subject_inds))
stc_alla2=range(len(subject_inds))
stc_allc2=range(len(subject_inds))

vertices_to = [np.arange(10242), np.arange(10242)]
b='gamma'
for ii, meg in enumerate(ll):
    print ii
    #fname = data_path + meg + 'mineMorphed_typlex_icaclean_word_Power_ratio_equalized_m500_700_' + b
    #for normori
    fname = data_path + meg + 'mineMorphed_ico_TypLex_ica_word_Power2_normori_ratio_equalized_m500_700_200hz_' +b
    stc_allc[ii] = mne.read_source_estimate(fname)
    stc_allc2[ii] = mne.read_source_estimate(fname)

    Matx0=np.zeros((20484,5))
    b2=0.05
    for cc in range(5):
        b1=b2+0.0001; b2=b1+0.1-0.0001
        Matx0[:,cc]=stc_allc[ii].copy().crop(b1,b2).mean().data.squeeze()
    matx_stcc = mne.SourceEstimate(Matx0, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
    stc_allc[ii]=matx_stcc-1
    #out_file0=data_path + meg + 'mineMorphed_typlex_icaclean_word_Power_ratio_equalized_meanresampled_50_550_100ms_0overlap_' +b
    #for normori
    out_file0=data_path + meg + 'mineMorphed_ico_TypLex_ica_word_Power2_normori_ratio_equalized_50_550_100ms_0ov_200hz_' +b
    matx_stcc.save(out_file0)
    
    fname = data_path + meg + 'mineMorphed_ico_TypLex_ica_nonword_Power2_normori_ratio_equalized_m500_700_200hz_' +b#'mineMorphed_typlex_icaclean_nonword_Power_ratio_equalized_m500_700_' + b  
    stc_alla[ii] = mne.read_source_estimate(fname)
    stc_alla2[ii] = mne.read_source_estimate(fname)

    Matx0=np.zeros((20484,5))
    b2=0.05
    for cc in range(5):
        b1=b2+0.0001; b2=b1+0.1-0.0001
        Matx0[:,cc]=stc_alla[ii].copy().crop(b1,b2).mean().data.squeeze()
    matx_stca = mne.SourceEstimate(Matx0, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
    stc_alla[ii]=matx_stca-1
    out_file0=data_path + meg + 'mineMorphed_ico_TypLex_ica_nonword_Power2_normori_ratio_equalized_50_550_100ms_0ov_200hz_' +b#'mineMorphed_typlex_icaclean_nonword_Power_ratio_equalized_meanresampled_50_550_100ms_0overlap_' +b
    matx_stca.save(out_file0)
    
    matx_stcs=np.subtract(matx_stcc,matx_stca)
    out_file0=data_path + meg + 'mineMorphed_ico_TypLex_ica_sub_Power2_normori_ratio_equalized_50_550_100ms_0ov_200hz_' +b#'mineMorphed_typlex_icaclean_subtract_Power_ratio_equalized_meanresampled_50_550_100ms_0overlap_' +b
    matx_stcs.save(out_file0)

stc_grandc=np.mean(stc_allc)
stc_grandc2=np.mean(stc_allc2)

#datag=np.log(stc_grand.data)
#datag_stc = mne.SourceEstimate(datag, vertices=stc_grand.vertno, tmin=stc_grand.tmin, tstep=stc_grand.tstep, subject='fsaverageect')
out_file=out_path + 'GrandAverage_mineMorphed_ico_TypLex_ica_word_Power2_normori_ratio_equalized_50_550_100ms_0ov_200hz_' +b#mineMorphed_typlex_icaclean_word_Power_ratio_equalized_meanresampled_50_550_100ms_0overlap_' +b
stc_grandc.save(out_file)

out_file=out_path + 'GrandAverage_mineMorphed_ico_TypLex_ica_word_Power2_normori_ratio_equalized_m500_700_200hz_' +b
stc_grandc2.save(out_file)

stc_granda=np.mean(stc_alla)
stc_granda2=np.mean(stc_alla2)

#datag=np.log(stc_grand.data)
#datag_stc = mne.SourceEstimate(datag, vertices=stc_grand.vertno, tmin=stc_grand.tmin, tstep=stc_grand.tstep, subject='fsaverageect')
out_file=out_path + 'GrandAverage_mineMorphed_ico_TypLex_ica_nonword_Power2_normori_ratio_equalized_50_550_100ms_0ov_200hz_' +b#mineMorphed_typlex_icaclean_nonword_Power_ratio_equalized_meanresampled_50_550_100ms_0overlap_' +b
stc_granda.save(out_file)
out_file=out_path + 'GrandAverage_mineMorphed_ico_TypLex_ica_nonword_Power2_normori_ratio_equalized_m500_700_200hz_' +b#mineMorphed_typlex_icaclean_nonword_Power_ratio_equalized_m500_700_' +b
stc_granda2.save(out_file)

stc_subtract=np.mean([stc_grandc,stc_granda])
out_file=out_path + 'GrandAverage_mineMorphed_ico_TypLex_ica_sub_Power2_normori_ratio_equalized_50_550_100ms_0ov_200hz_' +b#mineMorphed_typlex_icaclean_subtract_Power_ratio_equalized_meanresampled_50_550_100ms_0overlap_' +b 
stc_subtract.save(out_file)

stc_subtract2=np.mean([stc_grandc2,stc_granda2])
out_file=out_path + 'GrandAverage_mineMorphed_ico_TypLex_ica_sub_Power2_normori_ratio_equalized_m500_700_200hz_' +b#mineMorphed_typlex_icaclean_subtract_Power_ratio_equalized_m500_700_' +b 
stc_subtract2.save(out_file)

#fname=out_path + 'GrandAverage_preMorphed_SemLoc_icaclean_Concrete_Power_ratio_equalized_m500_700_theta'
#stcc=mne.read_source_estimate(fname)
#fname=out_path + 'GrandAverage_preMorphed_SemLoc_icaclean_Abstract_Power_ratio_equalized_m500_700_theta'
#stca=mne.read_source_estimate(fname)
#stcmt=np.mean([stcc,stca])
#
#fname=out_path + 'GrandAverage_preMorphed_SemLoc_icaclean_Concrete_Power_ratio_equalized_m500_700_alpha'
#stcc=mne.read_source_estimate(fname)
#fname=out_path + 'GrandAverage_preMorphed_SemLoc_icaclean_Abstract_Power_ratio_equalized_m500_700_alpha'
#stca=mne.read_source_estimate(fname)
#stcma=np.mean([stcc,stca])
#
#fname=out_path + 'GrandAverage_preMorphed_SemLoc_icaclean_Concrete_Power_ratio_equalized_m500_700_beta'
#stcc=mne.read_source_estimate(fname)
#fname=out_path + 'GrandAverage_preMorphed_SemLoc_icaclean_Abstract_Power_ratio_equalized_m500_700_beta'
#stca=mne.read_source_estimate(fname)
#stcmb=np.mean([stcc,stca])
#
#fname=out_path + 'GrandAverage_preMorphed_SemLoc_icaclean_Concrete_Power_ratio_equalized_m500_700_gamma'
#stcc=mne.read_source_estimate(fname)
#fname=out_path + 'GrandAverage_preMorphed_SemLoc_icaclean_Abstract_Power_ratio_equalized_m500_700_gamma'
#stca=mne.read_source_estimate(fname)
#stcmg=np.mean([stcc,stca])
#import matplotlib.pyplot as plt
#plt.figure(2)
#in_file1=out_path + 'GrandAverage_mineMorphed_ico_TypLex_ica_sub_Power2_normori_ratio_equalized_m500_700_200hz_theta'  
#stc1=mne.read_source_estimate(in_file1)
#in_file2=out_path + 'GrandAverage_mineMorphed_ico_TypLex_ica_sub_Power2_normori_ratio_equalized_m500_700_200hz_alpha' 
#stc2=mne.read_source_estimate(in_file2) 
#in_file3=out_path + 'GrandAverage_mineMorphed_ico_TypLex_ica_sub_Power2_normori_ratio_equalized_m500_700_200hz_beta'  
#stc3=mne.read_source_estimate(in_file3)
#in_file4=out_path + 'GrandAverage_mineMorphed_ico_TypLex_ica_sub_Power2_normori_ratio_equalized_m500_700_200hz_gamma'  
#stc4=mne.read_source_estimate(in_file4)
#plt.plot(stc1.times[50:220], stc1.data[:,50:220].mean(axis=0), label='Theta')
#plt.plot(stc2.times[50:220], stc2.data[:,50:220].mean(axis=0), label='Alpha')
#plt.plot(stc3.times[50:220], stc3.data[:,50:220].mean(axis=0), label='Beta')
#plt.plot(stc4.times[50:220], stc4.data[:,50:220].mean(axis=0), label='Gamma')
#plt.xlabel('Time (ms)')
#plt.ylabel('Power Ratio')
#plt.legend()
#plt.title('Average Source Evoked+Induced Power')
#plt.show()
#GrandAverage_preMorphed_typlex_icaclean_nonword_Power_ratio_equalized_m500_700_gamma-lh.stc

#
