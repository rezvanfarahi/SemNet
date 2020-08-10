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
sys.path.insert(1,'/imaging/rf02/scikit-learn-0.16.0')
sys.path.insert(1,'/home/rf02/.local/lib/python2.7/site-packages')

# for qsub
# (add ! here if needed) /imaging/local/software/anaconda/latest/x86_64/bin/python
#sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
###

import numpy as np
import mne

from mne.minimum_norm import ( make_inverse_operator, write_inverse_operator,read_inverse_operator)
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
   raw_fname = directory + 'typlex_pw_ssstf_fft_1_48_clean_ica_raw.fif'
   raw = mne.io.Raw(raw_fname, preload=True)

   print "Reading events"   
   events_fname = orig_path + meg + 'typlex_pw_raw_ssstnew.txt'
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
   if index < 3:
	    	word_ids=np.array([2,1]); print events[events[:,2]==1,2].shape+events[events[:,2]==2,2].shape
	    	nonword_ids=np.array([6,5]); print events[events[:,2]==5,2].shape+events[events[:,2]==6,2].shape
	    	events=mne.merge_events(events,word_ids,601,replace_events=True)
	    	events=mne.merge_events(events,nonword_ids,602,replace_events=True)
	    	event_id = {'Words': 601, 'PWords': 602}

   else:
	    	word_ids=np.array([8,7]); print events[events[:,2]==7,2].shape+events[events[:,2]==8,2].shape
	    	nonword_ids=np.array([6,5]); print events[events[:,2]==5,2].shape+events[events[:,2]==6,2].shape
	    	events=mne.merge_events(events,word_ids,601,replace_events=True)
	    	events=mne.merge_events(events,nonword_ids,602,replace_events=True)
	    	event_id = {'Words': 601, 'PWords': 602}
 
   exclude = []#raw.info['bads']
   include = []
   picks = mne.pick_types(raw.info, meg=True, eeg=True, eog=True, stim=False, include=include, exclude=exclude)
#   picks_eeg = mne.pick_types(raw.info, meg=False, eeg=True, eog=True, stim=False, include=include, exclude=exclude)
#   picks_meg = mne.pick_types(raw.info, meg=True, eeg=False, eog=True, stim=False, include=include, exclude=exclude)
   
   print "Epoching"
   if index in np.array([3,4,13]):
   	reject = dict(eeg=150e-6, grad=200e-12, mag=4e-12)# eog=150e-6,
        reject_eeg = dict(eeg=150e-6)
   else:
   	reject = dict(eeg=120e-6, grad=200e-12, mag=4e-12)# eog=150e-6,
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

   #read forward solution
   fname_fwd = data_path + subjects[index] + 'ico5_forward_5-3L-EMEG-fwd.fif'

   print "reading forward solution"
   forward_meeg = mne.read_forward_solution(fname_fwd, surf_ori=True)
#   forward_meg = mne.pick_types_forward(forward_meeg, meg=True, eeg=False)
   # Alternatively, you can just load a forward solution that is restricted
#   forward_eeg = mne.pick_types_forward(forward_meeg, meg=False, eeg=True)
   print "making noise covariance matrix"
   


   #Compute the noise covariance on baseline
   method_params=dict(diagonal_fixed=dict(grad=0.1, mag=0.1, eeg=0.1))
   noise_cov = mne.compute_covariance(epochs, method=['diagonal_fixed'], tmin=None, tmax=0., method_params=method_params, projs=raw.info['projs'])
#   noise_cov_meg = mne.compute_covariance(epochs_meg, method=['empirical','shrunk', 'diagonal_fixed'], tmin=None, tmax=0.)
#   noise_cov_eeg = mne.compute_covariance(epochs_eeg, method=['empirical','shrunk', 'diagonal_fixed'], tmin=None, tmax=0.)
   #regularize noise cov...
   #####noise_cov = mne.cov.regularize(noise_cov, epochs.info, mag=0.1, grad=0.1, eeg=0.1, proj=True) ##I checked with and without proj, both exactly the same in this case
   #mne.viz.plot_cov(noise_cov, raw.info, colorbar=True, proj=True)
  # make an M/EEG, MEG-only, and EEG-only inverse operators
   ####noise_cov_eeg = mne.cov.regularize(noise_cov_eeg, epochs_eeg.info, eeg=0.1, proj=True)
   ####noise_cov_meg = mne.cov.regularize(noise_cov_meg, epochs_meg.info, mag=0.1, grad=0.1, proj=True)
   info = evoked.info
   #noise_cov = mne.compute_covariance(epochs, tmin=None, tmax=0.)
#   info_eeg = evoked_eeg.info
#   info_meg = evoked_meg.info
   print "a"
   inverse_operator_meeg = make_inverse_operator(info, forward_meeg, noise_cov,loose=0.2, depth=None)
   #noise_cov = mne.cov.regularize(noise_cov, epochs.info, mag=0.1, grad=0.1, eeg=0.1, proj=True)
#   print "b"
#   inverse_operator_meg = make_inverse_operator(info_meg, forward_meg, noise_cov_meg, loose=0.2, depth=None)
#   print "c"
#   inverse_operator_eeg = make_inverse_operator(info_eeg, forward_eeg, noise_cov_eeg, loose=0.2, depth=None)
   fname_inv = data_path + meg + 'InvOp_ico5newreg_fft_1_48_clean_ica_EMEG-inv.fif'
   inverse_operator_old = read_inverse_operator(fname_inv)
   sub_src=inverse_operator_old['src']
#   #mne.add_source_space_distances(sub_src)
#   inverse_operator_meeg['src']=sub_src
#   sub_src=inverse_operator_meeg['src']
#   mne.add_source_space_distances(sub_src)
   inverse_operator_meeg['src']=sub_src
   #Write the inverse operator
   write_inverse_operator(data_path + subjects[index] + 'typlex_InvOp_ico5diagreg_fft_1_48_clean_ica_EMEG-inv.fif', inverse_operator_meeg)
#   write_inverse_operator(data_path + subjects[index] + 'InverseOperator_fft_1_48_clean_ica_MEG-inv.fif', inverse_operator_meg)
#   write_inverse_operator(data_path + subjects[index] + 'InverseOperator_fft_1_48_clean_ica_EEG-inv.fif', inverse_operator_eeg)
 

   ###end loop across subjects
