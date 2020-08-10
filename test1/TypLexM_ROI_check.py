# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 03:52:45 2015

@author: rf02
"""
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
import scipy


###############################################################################
out_path = '/imaging/rf02/TypLexMEG/icaanalysis_results/stc/GrandAverage/power/' # root directory for your MEG data
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
reject = dict(eeg=120e-6, eog=150e-6, grad=200e-12, mag=4e-12)
event_ids = {'cncrt_wrd': 1, 'abs_wrd': 2}
# frequency bands with variable number of cycles for wavelets
stim_delay = 0.034 # delay in s


tmin1=0
tstep1=50
stc_allc=range(17)
vertices_to = [np.arange(10242), np.arange(10242)]
stc_alla=range(17)
stc_alla2=range(17)
stc_allc2=range(17)
sdata=np.zeros((1,17))
cdatam=np.zeros((1,17))
adatam=np.zeros((1,17))
vertices_to = [np.arange(10242), np.arange(10242)]
for ii, meg in enumerate(ll):
    print ii
    fname = data_path + meg + 'firstMorphed_ico_typlex_ica_word_Source_Evoked_m500_700'  
    stcc = mne.read_source_estimate(fname)
    fname = data_path + meg + 'firstMorphed_ico_typlex_ica_nonword_Source_Evoked_m500_700'  
    stca = mne.read_source_estimate(fname)
    fname_label='/imaging/rf02/TypLexMEG/fsaverage/label/lh.WFA_DCM-lh.label'
    #fname_label='/imaging/rf02/TypLexMEG/createdlabels_SL/atlleft-lh.label'
    label1 = mne.read_label(fname_label)
    label1.values.fill(1.0)
    
    #label1=label1.morph(subject_from='fsaverage', subject_to='fsaverage', subjects_dir=data_path) 
    #stcc.resample(200); stca.resample(200)   
    cdata=stcc.in_label(label1).crop(.12,.150).data
    adata=stca.in_label(label1).crop(.12,.150).data
    cdatam[0,ii]=np.abs(cdata).max(axis=0).mean()
    adatam[0,ii]=np.abs(adata).max(axis=0).mean()
    sdata[0,ii]=cdatam[0,ii]-adatam[0,ii]

MT= mne.stats.ttest_1samp_no_p(sdata[0,:], sigma=1e-3, method='relative')
pval=scipy.stats.t.sf(np.abs(MT),16)*2
print pval
#path_con=data_path+'ROI_check_con.mat'
#scipy.io.savemat(path_con,{'cdatam':cdatam})
#path_abs=data_path+'ROI_check_abs.mat'
#scipy.io.savemat(path_abs,{'adatam':adatam})