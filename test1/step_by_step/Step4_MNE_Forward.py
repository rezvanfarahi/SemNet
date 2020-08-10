
"""
=========================================================
Noise covariance matrix and inverse operator
=========================================================


"""
# Authors: Alexandre Gramfort <gramfort@nmr.mgh.harvard.edu>
#
# License: BSD (3-clause)

print __doc__

import sys
sys.path.insert(1,'/imaging/local/software/mne_python/latest_v0.9full')
sys.path.insert(1,'/imaging/rf02/scikit-learn-0.15.0')
sys.path.insert(1,'/home/rf02/.local/lib/python2.7/site-packages')

# for qsub
# (add ! here if needed) /imaging/local/software/anaconda/latest/x86_64/bin/python
#sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
###

import numpy as np
import mne

from mne.minimum_norm import ( make_inverse_operator, write_inverse_operator)
import os


###############################################################################
data_path = '/imaging/rf02/TypLexMEG/' # root directory for your MEG data
orig_path='/group/erp/data/olaf.hauk/MEG/TypLexM/data/MF22_new/'
subjects_dir = '/imaging/rf02/TypLexMEG/'    # where your MRI subdirectories are
os.chdir('/imaging/rf02/TypLexMEG/')

# get indices for subjects to be processed from command line input
# 
print sys.argv
subject_inds = []#
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

subjects = []
for ss in subject_inds:
 subjects.append(list_all[ss])

print "subjects:"
print subjects



#subjects=[subjects[0]]
for index, meg in enumerate(subjects):
    print index #################################################################################
    #read concatenated events and raw file from disc
    directory = data_path + subjects[index]      
    print "reading raw file"
    raw_fname = directory + 'semloc_ssstf_fft_1_48_clean_ica_raw.fif'
    raw = mne.io.Raw(raw_fname, preload=True)
    print "Reading events"   
    events_fname = orig_path + meg + 'semloc_raw_ssstnew.txt'
    events = mne.read_events(events_fname)

    print "Reading evevokeds"   
    #fname_evoked = directory + 'semloc_ssstf_fft_1_48_clean_ica_raw-ave.fif'
    #evoked = mne.read_evokeds(fname_evoked, condition=0, baseline=(None, 0))
    stim_delay = 0.034 # delay in s
    events[:,0] += np.round( raw.info['sfreq']*stim_delay )
    ##################################################################################
    #extracting epochs
    tmin, tmax = -0.5, 0.7
    
    #Cue = {'EasyMul': 101, 'ComplexMul': 102, 'EasySeq': 103, 'ComplexSeq': 104}
    #Task = {'EasyMul': 131, 'ComplexMul': 132, 'EasySeq': 133, 'ComplexSeq': 134}
    event_id={'cncrt_wrd': 1, 'abs_wrd': 2}
    exclude = raw.info['bads']
    include = []
    picks = mne.pick_types(raw.info, meg=True, eeg=True, eog=True, stim=False, include=include, exclude=exclude)
    #   picks_eeg = mne.pick_types(raw.info, meg=False, eeg=True, eog=True, stim=False, include=include, exclude=exclude)
    #   picks_meg = mne.pick_types(raw.info, meg=True, eeg=False, eog=True, stim=False, include=include, exclude=exclude)
    print "Epoching"
    if index==3:
        reject = dict(eeg=150e-6, grad=200e-12, mag=4e-12)# eog=150e-6,
        reject_eeg = dict(eeg=150e-6)
    else:
        reject = dict(eeg=120e-6, grad=200e-12, mag=4e-12)# eog=150e-6
        reject_eeg = dict(eeg=120e-6)
    reject_meg = dict(grad=200e-12, mag=4e-12)   
    epochs = mne.Epochs(raw, events, event_id, tmin, tmax, picks=picks, proj=True, baseline=(None, 0), reject=reject)
    evoked=epochs.average()
#   epochs_eeg = mne.Epochs(raw, events, event_id, tmin, tmax, picks=picks_eeg, proj=True, baseline=(None, 0), reject=reject_eeg)
#   evoked_eeg=epochs_eeg.average()
#   epochs_meg = mne.Epochs(raw, events, event_id, tmin, tmax, picks=picks_meg, proj=True, baseline=(None, 0), reject=reject_meg)
#   evoked_meg=epochs_meg.average()
   #epochs=mne.epochs.equalize_epoch_counts(epochs)
   ##################################################################################

   #do forward solution
    fname_fwd = data_path + subjects[index] + 'mnepython_forward_5-3L-EMEG-fwd.fif'
    bem_file=directory+'/bem/'+subjects[index]+'-5120-5120-5120-bem-sol.fif'
    src_file=directory+'/bem/'+subjects[index]+'-5-src.fif'
    mri_file=directory+'/T1-neuromag/sets/COR.fif'
    print "forward solution"
    forward_meeg = mne.do_forward_solution(subject=subjects[index],meas=evoked,fname=fname_fwd,src=src_file,spacing='ico5',mindist=5.0,bem=bem_file,mri=mri_file, eeg=True, meg=True, fixed=False, overwrite=True, subjects_dir=data_path)
    