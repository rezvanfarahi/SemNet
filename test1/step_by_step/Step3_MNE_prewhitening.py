# -*- coding: utf-8 -*-
"""
Created on Mon May 11 13:14:33 2015

@author: rf02
"""
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
sys.path.insert(1,'/imaging/local/software/mne_python/latest_v0.9')
sys.path.insert(1,'/imaging/rf02/scikit-learn-0.15.0')
sys.path.insert(1,'/home/rf02/.local/lib/python2.7/site-packages')
import sklearn
reload(sklearn)

# for qsub
# (add ! here if needed) /imaging/local/software/anaconda/latest/x86_64/bin/python
#sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
###

import numpy as np
import mne
from mne.io import Raw
import os
import scipy.io
from mne.cov import _get_whitener_data as mne_gwd
from mne.cov import compute_whitener
import time

###############################################################################
ts=time.time()
data_path = '/imaging/rf02/TypLexMEG/' # root directory for your MEG data
orig_path='/group/erp/data/olaf.hauk/MEG/TypLexM/data/MF22_new/'
subjects_dir = '/imaging/rf02/TypLexMEG/'    # where your MRI subdirectories are
os.chdir('/imaging/rf02/TypLexMEG/')

# get indices for subjects to be processed from command line input
# 
print sys.argv
subject_inds = [0]#
for ss in sys.argv[1:]:
    subject_inds.append( int( ss ) )

#subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
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
 
   exclude = []#raw.info['bads']
   include = []
   picks = mne.pick_types(raw.info, meg=True, eeg=True, eog=False, stim=False)
   
   print "Epoching"
   if index==3:
   	reject = dict(eeg=150e-6, grad=200e-12, mag=4e-12)# eog=150e-6,
   else:
   	reject = dict(eeg=120e-6, grad=200e-12, mag=4e-12)# eog=150e-6, 
   epochs = mne.Epochs(raw, events, event_id, tmin, tmax, picks=picks, proj=True, baseline=(None, 0), reject=reject)
   evoked=epochs.average()
   #epochs=mne.epochs.equalize_epoch_counts(epochs)
   ##################################################################################
#
   #read forward solution
   fname_fwd = data_path + subjects[index] + 'forward_5-3L-EMEG-fwd.fif'

   print "reading forward solution"
   forward_meeg = mne.read_forward_solution(fname_fwd, surf_ori=True)
   ec=epochs['cncrt_wrd']; ea=epochs['abs_wrd']

   #Compute the noise covariance on baseline

   noise_cov = mne.compute_covariance(epochs, method=['empirical','shrunk', 'diagonal_fixed'], tmin=None, tmax=0.)
   #noise_cov_a = mne.compute_covariance(ea, method=['empirical','shrunk', 'diagonal_fixed'], tmin=None, tmax=0.)

   info = epochs.info
   #ch_names = [info['ch_names'][k] for k in picks]
   #noise_cov = pick_channels_cov(noise_cov, include=ch_names, exclude=[])

   scalings_ = dict(mag=1e15, grad=1e13, eeg=1e6)

   W = compute_whitener(noise_cov, info, rank=None, scalings=scalings_)[0]
   #Wa = compute_whitener(noise_cov_a, info, rank=None, scalings=scalings_)[0]
   ecd=ec.get_data(); ead=ea.get_data()
   #ed=epochs.get_data()
   ecdw=ecd.copy(); eadw=ead.copy()
   for cntt in range(ecd.shape[0]):
       ecdw[cntt,:,:]=np.dot(W, ecd[cntt,:,:])
   for cntt in range(ead.shape[0]):
       eadw[cntt,:,:]=np.dot(W, ead[cntt,:,:])
   mne.EpochsArray(ecdw, ec.info, tmin, tmax, event_id['cncrt_wrd'])
       
    
#   ecdw=np.dot(Wc, ecd); eadw=np.dot(Wa, eca)
#   epochs_concrete=mne.EpochsArray(ecdw, ec.info, events)
#   raw[picks]=raw_dat
#   data_out_path=directory + 'semloc_ssstf_fft_1_48_clean_ica_whitened_raw.fif'
#   raw_new.save(data_out_path, overwrite=True)
#   #M=mne_gwd(evoked.info, noise_cov, picks, diag=False, rank=None, scalings=None, verbose=None)
   te=time.time()-ts
   print te
   
#   def epochs_whiten(epochs, W=None, verbose=None):
#        
#        for cntt in range(epochs._data.shape[0]):
#            epochs._data[cntt,:,:]=np.dot(W, epochs._data[cntt,:,:])
#        return epochs