"""
=========================================================================
Compute source space connectivity and visualize it using a circular graph
=========================================================================

This example computes the all-to-all connectivity between 68 regions in
source space based on dSPM inverse solutions and a FreeSurfer cortical
parcellation. The connectivity is visualized using a circular graph which
is ordered based on the locations of the regions.
"""

# Authors: Martin Luessi <mluessi@nmr.mgh.harvard.edu>
#          Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
#          Nicolas P. Rougier (graph code borrowed from his matplotlib gallery)
#
# License: BSD (3-clause)

print(__doc__)
import sys

sys.path.append('/imaging/local/software/python_packages/nibabel/1.3.0')
sys.path.append('/imaging/local/software/python_packages/pysurfer/v0.4')
#sys.path.append('/imaging/local/software/python_packages/scikit-learn/v0.14.1')
#sys.path.append('/imaging/local/software/mne_python/latest')

# for qsub
# (add ! here if needed) 
sys.path.append('/imaging/local/software/EPD/latest/x86_64/bin/python')
sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
###

import os
import numpy as np
import mne
from mne import io
from mne.io import Raw
import operator
from mne.minimum_norm import (read_inverse_operator, make_inverse_operator, write_inverse_operator, apply_inverse, apply_inverse_epochs)
from mne.minimum_norm.inverse import (prepare_inverse_operator)
import sklearn
import scipy.io
from mne import find_events
from mne.connectivity import seed_target_indices, spectral_connectivity
import numpy as np
import mne
from mne.datasets import sample
from mne.io import Raw
from mne.minimum_norm import apply_inverse_epochs, read_inverse_operator
from mne.connectivity import spectral_connectivity
from mne.viz import circular_layout, plot_connectivity_circle


data_path = '/imaging/rf02/TypLexMEG/' # root directory for your MEG data
subjects_dir = '/home/rf02/rezvan/TypLexMEG/'    # where your MRI subdirectories are

event_path = '/imaging/rf02/TypLexMEG/'    # where event files are

label_path = '/home/rf02/rezvan/TypLexMEG/createdlabels_SL/'

inv_path = '/imaging/rf02/TypLexMEG/'
inv_fname = 'InverseOperator_EMEG-inv.fif'

# get indices for subjects to be processed from command line input
# 
print sys.argv
subject_inds = []
for ss in sys.argv[1:]:
   subject_inds.append( int( ss ) )

subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
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

# subjects names used for MRI data
subjects=['SS1_Meg10_0378',
'SS2_Meg10_0390',
'SS3_Meg11_0026',
'SS4_Meg11_0050',
'SS5_Meg11_0052',
#'SS6_Meg11_0069',
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
Dmat=np.zeros((16,4,4512))
for cnt1, meg in enumerate(ll):
	print meg
	
	fname_inv = inv_path + meg + 'TypLex_PW_InverseOperator_EMEG-inv.fif'
	fname_raw = data_path + meg + 'typlex_pw_ssstf_fft49_raw.fif'
	fname_event = event_path + meg + 'typlex_pw_ssstf_fft49_raw-eve.fif'

	# Load data
	inverse_operator = read_inverse_operator(fname_inv)
	raw = Raw(fname_raw)
	events = mne.read_events(fname_event)

	# Add a bad channel
	#raw.info['bads'] += ['MEG 2443']

	# Pick MEG channels
	picks = mne.pick_types(raw.info, meg=True, eeg=True, stim=False, eog=True,
		               exclude='bads')

	# Define epochs for left-auditory condition
	tmin, tmax = -0.5, 0.7
	reject = dict(eeg=120e-6, eog=150e-6, grad=200e-12, mag=4e-12)
	if subject_inds[cnt1] < 3:
    		event_ids = {'TypWords': 2, 'AtypWords': 1, 'TypPWords': 6, 'AtypPWords': 5}
    	else:
    		event_ids = {'TypWords': 8, 'AtypWords': 7, 'TypPWords': 6, 'AtypPWords': 5}

	# frequency bands with variable number of cycles for wavelets
	stim_delay = 0.034 # delay in s
	
	epochs = mne.Epochs(raw, events, event_ids, tmin, tmax, picks=picks,
		            baseline=(None, 0), reject=reject)
	epochs_TypWords = epochs['TypWords']
	epochs_AtypWords = epochs['AtypWords']
	epochs_TypPWords = epochs['TypPWords']
	epochs_AtypPWords = epochs['AtypPWords']

	# Compute inverse solution and for each epoch. By using "return_generator=True"
	# stcs will be a generator object instead of a list.
	snr = 1.0  # use lower SNR for single epochs
	lambda2 = 1.0 / snr ** 2
	method = "MNE"  # use dSPM method (could also be MNE or sLORETA)
	stcs_TypWords = apply_inverse_epochs(epochs_TypWords, inverse_operator, lambda2, method,
		                    pick_ori="normal", return_generator=True)
	stcs_AtypWords = apply_inverse_epochs(epochs_AtypWords, inverse_operator, lambda2, method,
		                    pick_ori="normal", return_generator=True)
	stcs_TypPWords = apply_inverse_epochs(epochs_TypPWords, inverse_operator, lambda2, method,
		                    pick_ori="normal", return_generator=True)
	stcs_AtypPWords = apply_inverse_epochs(epochs_AtypPWords, inverse_operator, lambda2, method,
		                    pick_ori="normal", return_generator=True)

	# Get labels for FreeSurfer 'aparc' cortical parcellation with 34 labels/hemi
	labels = mne.read_labels_from_annot(subject=subjects[subject_inds[cnt1]], parc='aparc',
		                            subjects_dir=subjects_dir)

	b1=[0,1,12,13,14,15,16,17,18,19,20,21,22,23,28,29,30,31,32,33,34,35,36,37,38,39, 40,41,42,43,46,47,48,49,50,51,54,55,58,59,62,63,64,65,66,67,68,69]
	labels=[labels[k1] for k1 in b1]

	label_colors = [label.color for label in labels]

	# Average the source estimates within each label using sign-flips to reduce
	# signal cancellations, also here we return a generator
	src = inverse_operator['src']
	label_ts_TypWords = mne.extract_label_time_course(stcs_TypWords, labels, src, mode='mean_flip',
		                                 return_generator=True)
	label_ts_AtypWords = mne.extract_label_time_course(stcs_AtypWords, labels, src, mode='mean_flip',
		                                 return_generator=True)
	label_ts_TypPWords = mne.extract_label_time_course(stcs_TypPWords, labels, src, mode='mean_flip',
		                                 return_generator=True)
	label_ts_AtypPWords = mne.extract_label_time_course(stcs_AtypPWords, labels, src, mode='mean_flip',
		                                 return_generator=True)

	# Now we are ready to compute the connectivity in the alpha band. Notice
	# from the status messages, how mne-python: 1) reads an epoch from the raw
	# file, 2) applies SSP and baseline correction, 3) computes the inverse to
	# obtain a source estimate, 4) averages the source estimate to obtain a
	# time series for each label, 5) includes the label time series in the
	# connectivity computation, and then moves to the next epoch. This
	# behaviour is because we are using generators and allows us to
	# compute connectivity in computationally efficient manner where the amount
	# of memory (RAM) needed is independent from the number of epochs.
	fmin = (4.)#,8., 13., 31.)
	fmax = (8.)#, 12., 30., 45.)
	tmin2 = 0.15
	tmax2 = 0.45
	sfreq = raw.info['sfreq']  # the sampling frequency
	con_methods = ['coh']
	con_TypWords, freqs_TypWords, times_TypWords, n_epochs_TypWords, n_tapers_TypWords = spectral_connectivity(label_ts_TypWords, method=con_methods, mode='multitaper', sfreq=sfreq, fmin=fmin, fmax=fmax, faverage=True, tmin=tmin2, tmax=tmax2, mt_adaptive=True, n_jobs=2)
	con_AtypWords, freqs_AtypWords, times_AtypWords, n_epochs_AtypWords, n_tapers_AtypWords = spectral_connectivity(label_ts_AtypWords,method=con_methods, mode='multitaper', sfreq=sfreq, fmin=fmin, fmax=fmax, faverage=True, tmin=tmin2, tmax=tmax2, mt_adaptive=True, n_jobs=2)

	con_TypPWords, freqs_TypPWords, times_TypPWords, n_epochs_TypPWords, n_tapers_TypPWords = spectral_connectivity(label_ts_TypPWords, method=con_methods, mode='multitaper', sfreq=sfreq, fmin=fmin, fmax=fmax, faverage=True, tmin=tmin2, tmax=tmax2, mt_adaptive=True, n_jobs=2)
	con_AtypPWords, freqs_AtypPWords, times_AtypPWords, n_epochs_AtypPWords, n_tapers_AtypPWords = spectral_connectivity(label_ts_AtypPWords,method=con_methods, mode='multitaper', sfreq=sfreq, fmin=fmin, fmax=fmax, faverage=True, tmin=tmin2, tmax=tmax2, mt_adaptive=True, n_jobs=2)

	# con is a 3D array, get the connectivity for the first (and only) freq. band
	# for each method
	"""
	con_res_TypWords = dict()
	for method, c in zip(con_methods, con_TypWords):
	    con_res_TypWords[method] = c[:, :, 0]

	con_res_AtypWords = dict()
	for method, c in zip(con_methods, con_AtypWords):
	    con_res_AtypWords[method] = c[:, :, 0]
	
	cc_TypWords=con_res_TypWords['plv']
	cc_AtypWords=con_res_AtypWords['plv']
	"""
	"""
	cc_TypWords=con_TypWords
	cc_AtypWords=con_AtypWords
	#kk=0
	for ii in range(cc_TypWords.shape[0]-1):
		for jj in range(cc_TypWords.shape[1]-1):
			if ii<jj:
				#kk=kk+1
				cc_TypWords[ii,jj,:]=cc_TypWords[jj,ii,:]

	TypWords_lh=cc_TypWords[::2,::2,:]
	TypWords_lhmean=TypWords_lh.sum(axis=0).sum(axis=0)/((TypWords_lh.shape[0]-1)*(TypWords_lh.shape[1]-1))
	TypWords_rh=cc_TypWords[1::2,1::2]
	TypWords_rhmean=TypWords_rh.sum(axis=0).sum(axis=0)/((TypWords_rh.shape[0]-1)*(TypWords_rh.shape[1]-1))
	
	cc2=[]
	#kk=0
	for ii in range(cc_AtypWords.shape[0]-1):
		for jj in range(cc_AtypWords.shape[1]-1):
			if ii<jj:
				#kk=kk+1
				cc_AtypWords[ii,jj]=cc_AtypWords[jj,ii]

	AtypWords_lh=cc_AtypWords[::2,::2]
	AtypWords_lhmean=AtypWords_lh.sum(axis=0).sum(axis=0)/((AtypWords_lh.shape[0]-1)*(AtypWords_lh.shape[1]-1))
	AtypWords_rh=cc_AtypWords[1::2,1::2]
	AtypWords_rhmean=AtypWords_rh.sum(axis=0).sum(axis=0)/((AtypWords_rh.shape[0]-1)*(AtypWords_rh.shape[1]-1))
	"""
	cc_TypWords=con_TypWords[con_TypWords!=0]
	cc_AtypWords=con_AtypWords[con_AtypWords!=0]

	cc_TypPWords=con_TypPWords[con_TypPWords!=0]
	cc_AtypPWords=con_AtypPWords[con_AtypPWords!=0]

	Dmat[cnt1,:,:]=[cc_TypWords,cc_AtypWords,cc_TypPWords,cc_AtypPWords]


factor_levels = [2, 2]  # number of levels in each factor
effects = 'A*B' 
fvals, pvals = mne.stats.f_twoway_rm(Dmat, factor_levels, effects=effects)
reject, pval_corrected=mne.stats.fdr_correction(pvals[1,:], alpha=0.05, method='indep')
#aa0=pval_corrected[0,reject[0,:]]
		
	
