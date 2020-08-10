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
sys.path.append('/imaging/local/software/python_packages/nibabel/1.3.0')
sys.path.append('/imaging/local/software/python_packages/pysurfer/v0.4')
sys.path.append('/imaging/local/software/python_packages/scikit-learn/v0.14.1')
sys.path.append('/imaging/local/software/mne_python/latest')

# for qsub
# (add ! here if needed) /imaging/local/software/anaconda/latest/x86_64/bin/python
sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
###

import matplotlib.pyplot as plt
import numpy as np
import mne
from mne import io
from mne.io import Raw

from mne import (fiff, read_evokeds)
from mne.preprocessing.ica import ICA
import scipy.io as sio
import operator

from mne.minimum_norm import (read_inverse_operator, make_inverse_operator, write_inverse_operator, 
                              apply_inverse)
from mne.minimum_norm.inverse import (prepare_inverse_operator, _assemble_kernel)

from surfer import Brain
from surfer.io import read_stc
import logging
import os
import sklearn
import scipy.io
from mne import filter
from mne import find_events
from mne.epochs import combine_event_ids


###############################################################################
data_path = '/imaging/rf02/TypLexMEG/' # root directory for your MEG data
subjects_dir = '/imaging/rf02/TypLexMEG/'    # where your MRI subdirectories are
#os.chdir('/home/rf02/rezvan/test1')

# get indices for subjects to be processed from command line input
# 
print sys.argv
subject_inds=[]
for ss in sys.argv[1:]:
    subject_inds.append( int( ss ) )
subject_inds = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]

#subject_inds=[0,1]
print "subject_inds:"
print subject_inds
print "No rejection"
list_all = ['meg10_0378/101209/', 
'meg10_0390/101214/',
'meg11_0026/110223/', 
'meg11_0050/110307/', 
'meg11_0052/110307/', 
#'meg11_0069/110315/', 
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

ll = []
for ss in subject_inds:
 ll.append(list_all[ss])

print "subjects:"
print ll


for ii, meg in enumerate(ll):
   print ii, meg
   if subject_inds[ii] < 3:
    	event_ids = {'TypWords': 2, 'AtypWords': 1, 'TypPWords': 6, 'AtypPWords': 5}
   else:
    	event_ids = {'TypWords': 8, 'AtypWords': 7, 'TypPWords': 6, 'AtypPWords': 5}
   
   print "reading raw file"
   raw_fname = data_path + meg + 'typlex_pw_ssstf_fft49_raw.fif'
   raw = mne.io.Raw(raw_fname, preload=True)

   print "Reading events"   
   events_fname = data_path + meg + 'typlex_pw_ssstf_fft49_raw-eve.fif'
   events = mne.read_events(events_fname)

   print "Reading evevokeds"   
   fname_evoked = data_path + meg + 'typlex_pw_ssstf_fft49_raw-ave.fif'
   evoked = mne.read_evokeds(fname_evoked, condition=0, baseline=(None, 0))
   
   stim_delay = 0.034 # delay in s
   events[:,0] += np.round( raw.info['sfreq']*stim_delay )

   ##################################################################################
   
   #extracting epochs
   tmin, tmax = -0.5, 0.7

   #Cue = {'EasyMul': 101, 'ComplexMul': 102, 'EasySeq': 103, 'ComplexSeq': 104}
   #Task = {'EasyMul': 131, 'ComplexMul': 132, 'EasySeq': 133, 'ComplexSeq': 134}
 
   exclude = raw.info['bads']
   include = []
   picks = mne.pick_types(raw.info, meg=True, eeg=True, eog=True, stim=True, include=include, exclude=exclude)
   
   print "Epoching"
   reject = dict(eeg=120e-6, eog=150e-6, grad=200e-12, mag=4e-12)
   epochs = mne.Epochs(raw, events, event_ids, tmin, tmax, picks=picks, proj=True, baseline=(None, 0), reject=reject)
   ##################################################################################

   #read forward solution
   fname_fwd = data_path + meg + 'forward_5-3L-EMEG-fwd.fif'

   print "reading forward solution"
   forward_meeg = mne.read_forward_solution(fname_fwd, surf_ori=True)
   forward_meg = mne.pick_types_forward(forward_meeg, meg=True, eeg=False)
# Alternatively, you can just load a forward solution that is restricted
   forward_eeg = mne.pick_types_forward(forward_meeg, meg=False, eeg=True)


   #Compute the noise covariance on baseline
   noise_cov = mne.compute_covariance(epochs, tmin=None, tmax=0.)
   
   #regularize noise cov...
   noise_cov = mne.cov.regularize(noise_cov, epochs.info, mag=0.05, grad=0.05, eeg=0.1, proj=True)

  # make an M/EEG, MEG-only, and EEG-only inverse operators
   info = evoked.info
   print "a"
   inverse_operator_meeg = make_inverse_operator(info, forward_meeg, noise_cov,
		                              loose=0.2, depth=0.8)
   print "b"
   inverse_operator_meg = make_inverse_operator(info, forward_meg, noise_cov,
		                             loose=0.2, depth=0.8)
   print "c"
   inverse_operator_eeg = make_inverse_operator(info, forward_eeg, noise_cov,
		                             loose=0.2, depth=0.8)

   #Write the inverse operator
   write_inverse_operator(data_path + meg + 'TypLex_PW_InverseOperator_EMEG-inv.fif', inverse_operator_meeg)
   write_inverse_operator(data_path + meg + 'TypLex_PW_InverseOperator_MEG-inv.fif', inverse_operator_meg)
   write_inverse_operator(data_path + meg + 'TypLex_PW_InverseOperator_EEG-inv.fif', inverse_operator_eeg)
 

   ###end loop across subjects

