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


###############################################################################
out_path = '/imaging/rf02/TypLexMEG/icaanalysis_results/stc/GrandAverage/evoked/' # root directory for your MEG data
subjects_dir = '/imaging/rf02/TypLexMEG/'    # where your MRI subdirectories are
# where event-files are
data_path = '/imaging/rf02/TypLexMEG/'    # where event files are

label_path = '/home/rf02/rezvan/TypLexMEG/createdlabels_SL/'

inv_path = '/imaging/rf02/TypLexMEG/'
inv_fname = 'InverseOperator_oldreg_EMEG-inv.fif'

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
reject = dict(eeg=120e-6, eog=150e-6, grad=200e-12, mag=4e-12)
event_ids = {'cncrt_wrd': 1, 'abs_wrd': 2}
# frequency bands with variable number of cycles for wavelets
stim_delay = 0.034 # delay in s
frequencies = np.arange(8, 23, 1)
n_cycles = frequencies / float(7)
n_cycles[frequencies<=20] = 2
n_cycles[frequencies<=20] += np.arange(0,1,0.077)
    # n_cycles[np.logical_and((frequencies>10), (frequencies<=20))] = 3

tmin1=50
tstep1=100
stc_allc=range(17)
vertices_to = [np.arange(10242), np.arange(10242)]
stc_alla=range(17)
vertices_to = [np.arange(10242), np.arange(10242)]
for ii, meg in enumerate(ll):
    print ii
    fname = data_path + meg + 'firstMorphed_SemLoc_icaclean_oldreg_Concrete_Source_Evoked_m500_700'  
    stc_allc[ii] = mne.read_source_estimate(fname)
    Matx0=np.zeros((20484,5))
    b2=0.05
    for cc in range(5):
        b1=b2+0.0001; b2=b1+0.1-0.0001
        Matx0[:,cc]=stc_allc[ii].copy().crop(b1,b2).mean().data.squeeze()
    matx_stcc = mne.SourceEstimate(Matx0, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
    stc_allc[ii]=matx_stcc
    out_file0=data_path + meg + 'firstMorphed_SemLoc_icaclean_oldreg_Concrete_Source_Evoked_meanresampled_50_550_100ms_0overlap'
    matx_stcc.save(out_file0)
    
    fname = data_path + meg + 'firstMorphed_SemLoc_icaclean_oldreg_Abstract_Source_Evoked_m500_700'  
    stc_alla[ii] = mne.read_source_estimate(fname)
    Matx0=np.zeros((20484,5))
    b2=0.05
    for cc in range(5):
        b1=b2+0.0001; b2=b1+0.1-0.0001
        Matx0[:,cc]=stc_alla[ii].copy().crop(b1,b2).mean().data.squeeze()
    matx_stca = mne.SourceEstimate(Matx0, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
    stc_alla[ii]=matx_stca
    out_file0=data_path + meg + 'firstMorphed_SemLoc_icaclean_oldreg_Abstract_Source_Evoked_meanresampled_50_550_100ms_0overlap'
    matx_stca.save(out_file0)
    
    matx_stcs=np.subtract(matx_stcc,matx_stca)
    out_file0=data_path + meg + 'firstMorphed_SemLoc_icaclean_oldreg_Subtract_Source_Evoked_meanresampled_50_550_100ms_0overlap'
    matx_stcs.save(out_file0)

stc_grandc=np.mean(stc_allc)
#datag=np.log(stc_grand.data)
#datag_stc = mne.SourceEstimate(datag, vertices=stc_grand.vertno, tmin=stc_grand.tmin, tstep=stc_grand.tstep, subject='fsaverageect')
out_file=out_path + 'GrandAverage_firstmorphed_SemLoc_icaclean_oldreg_Concrete_Source_Evoked_meanresampled_50_550_100ms_0overlap'
stc_grandc.save(out_file)

stc_granda=np.mean(stc_alla)
#datag=np.log(stc_grand.data)
#datag_stc = mne.SourceEstimate(datag, vertices=stc_grand.vertno, tmin=stc_grand.tmin, tstep=stc_grand.tstep, subject='fsaverageect')
out_file=out_path + 'GrandAverage_firstmorphed_SemLoc_icaclean_oldreg_Abstract_Source_Evoked_meanresampled_50_550_100ms_0overlap'
stc_granda.save(out_file)


stc_subtract=np.subtract(stc_grandc,stc_granda)
out_file=out_path + 'GrandAverage_firstmorphed_SemLoc_icaclean_oldreg_Subtract_Source_Evoked_meanresampled_50_550_100ms_0overlap' 
stc_subtract.save(out_file)

"""
import matplotlib.pyplot as plt
in_file1=out_path + 'GrandAverage_SemLoc_icaclean_oldreg_Subtract_SourcePower_ratio_m500_700_theta'  
stc1=mne.read_source_estimate(in_file1)
in_file2=out_path + 'GrandAverage_SemLoc_icaclean_oldreg_Subtract_SourcePower_ratio_m500_700_alpha' 
stc2=mne.read_source_estimate(in_file2) 
in_file3=out_path + 'GrandAverage_SemLoc_icaclean_oldreg_Subtract_SourcePower_ratio_m500_700_beta'  
stc3=mne.read_source_estimate(in_file3)
in_file4=out_path + 'GrandAverage_SemLoc_icaclean_oldreg_Subtract_SourcePower_ratio_m500_700_gamma'  
stc4=mne.read_source_estimate(in_file4)
plt.plot(stc1.times[100:1100], stc1.data[:,100:1100].mean(axis=0), label='theta')
plt.plot(stc2.times[100:1100], stc2.data[:,100:1100].mean(axis=0), label='alpha')
plt.plot(stc3.times[100:1100], stc3.data[:,100:1100].mean(axis=0), label='beta')
plt.plot(stc4.times[100:1100], stc4.data[:,100:1100].mean(axis=0), label='gamma')
plt.xlabel('Time (ms)')
plt.ylabel('Power')
plt.legend()
plt.title('Mean source induced power Subtract')
plt.show()
"""


