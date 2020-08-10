# Author: Martin Luessi <mluessi@nmr.mgh.harvard.edu>
#
# License: BSD (3-clause)

print(__doc__)

import sys

sys.path.append('/imaging/local/software/python_packages/nibabel/1.3.0')
sys.path.append('/imaging/local/software/python_packages/pysurfer/v0.4')
#sys.path.append('/imaging/local/software/python_packages/scikit-learn/v0.14.1')
#sys.path.append('/imaging/local/software/mne_python/latest')

# for qsub
# (add ! here if needed) /imaging/local/software/anaconda/latest/x86_64/bin/python
sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
###
import os
#import matplotlib.pyplot as plt
import numpy as np
import mne
from mne import io
from mne.io import Raw
from mne import (fiff, equalize_channels)
from mne.preprocessing.ica import ICA
import scipy.io as sio
import operator

from mne.minimum_norm import (read_inverse_operator, make_inverse_operator, write_inverse_operator, apply_inverse, read_inverse_operator)
from mne.minimum_norm.inverse import (prepare_inverse_operator, _assemble_kernel)

from surfer.io import read_stc
import logging
import os
import sklearn
import scipy.io
from mne import filter
from mne import find_events
from mne.epochs import combine_event_ids
from mne.minimum_norm import (apply_inverse, apply_inverse_epochs,
                              read_inverse_operator)
from mne.connectivity import seed_target_indices, spectral_connectivity

data_path = '/home/rf02/rezvan/TypLexMEG/' # root directory for your MEG data
os.chdir(data_path)
event_path = '/home/rf02/rezvan/TypLexMEG/'    # where event files are
label_path = '/group/erp/data/olaf.hauk/MEG/TypLexM/labels/createdlabels_SL/'
inv_path = '/home/rf02/rezvan/TypLexMEG/'
subjects_dir = data_path 
print sys.argv
subject_inds = []
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




#label_name_lh = ['atlleft-lh', 'atlright-rh']


# Labels/ROIs
labellist = ['atlleft-lh', 'atlright-rh', 'medialtempright-rh','medialtempleft-lh']

ii=-1
for meg in ll:
	tmin, tmax = -0.2, 0.5
	reject = dict(eeg=120e-6, eog=150e-6, grad=200e-12, mag=4e-12)
	event_ids = {'cncrt_wrd': 1, 'abs_wrd': 2}
	# frequency bands with variable number of cycles for wavelets
	stim_delay = 0.034 # delay in s
	frequencies = np.arange(8, 45, 1)
	n_cycles = frequencies / float(7)
	n_cycles[frequencies<=20] = 2
	n_cycles[frequencies<=20] += np.arange(0,1,0.077)
	    # n_cycles[np.logical_and((frequencies>10), (frequencies<=20))] = 3

	method = "MNE"  # use dSPM method (could also be MNE or sLORETA)

	ii=ii+1
	fname_inv = inv_path + meg + 'InverseOperator_EMEG-inv.fif'
	fname_raw = data_path + meg + '/semloc_raw_ssstf_raw.fif'
	fname_event = event_path + meg + 'semloc_raw_ssstf-eve.fif'
# Load data
	inverse_operator = read_inverse_operator(fname_inv)
	raw = Raw(fname_raw)
	events = mne.read_events(fname_event)

		# pick MEG channels
	picks = mne.pick_types(raw.info, meg=True, eeg=True, stim=False, eog=True, exclude='bads')

		# Read epochs
	epochs = mne.Epochs(raw, events, event_ids, tmin, tmax, picks=picks, baseline=(None, 0), reject=reject)
	epochs = epochs['cncrt_wrd']

	# First, we find the most active vertex in the left auditory cortex, which
	# we will later use as seed for the connectivity computation
	snr = 3.0
	lambda2 = 1.0 / snr ** 2
	evoked = epochs.average()
	

	for label_name in labellist:
		print label_name
		tmin, tmax = -0.2, 0.5
		# frequency bands with variable number of cycles for wavelets
		stim_delay = 0.034 # delay in s
		frequencies = np.arange(8, 45, 1)
		n_cycles = frequencies / float(7)
		n_cycles[frequencies<=20] = 2
		n_cycles[frequencies<=20] += np.arange(0,1,0.077)
		    # n_cycles[np.logical_and((frequencies>10), (frequencies<=20))] = 3

		method = "MNE"  # use dSPM method (could also be MNE or sLORETA)

		fname_label = label_path + subjects[ii] + '/' + label_name + '.label'
		label1 = mne.read_label(fname_label)
		
# Restrict the source estimate to the label in the left auditory cortex
		stc = apply_inverse(evoked, inverse_operator, lambda2, method,
			    pick_ori="normal")
		stc_label = stc.in_label(label1)

# Find number and index of vertex with most power
		src_pow = np.sum(stc_label.data ** 2, axis=1)
		if label_name[-2] is 'l':
			print label_name
			seed_vertno = stc_label.vertno[0][np.argmax(src_pow)]
			seed_idx = np.searchsorted(stc.vertno[0], seed_vertno)  # index in original stc

		else:
			print label_name
			seed_vertno = stc_label.vertno[1][np.argmax(src_pow)]
			seed_idx = np.searchsorted(stc.vertno[1], seed_vertno)  # index in original stc

# Generate index parameter for seed-based connectivity analysis
		n_sources = stc.data.shape[0]
		indices = seed_target_indices([seed_idx], np.arange(n_sources))

# Compute inverse solution and for each epoch. By using "return_generator=True"
		# stcs will be a generator object instead of a list. This allows us so to
		# compute the coherence without having to keep all source estimates in memory.

		snr = 1.0  # use lower SNR for single epochs
		lambda2 = 1.0 / snr ** 2
		stcs = apply_inverse_epochs(epochs, inverse_operator, lambda2, method, pick_ori="normal", return_generator=True)

# Now we are ready to compute the coherence in the alpha and beta band.
# fmin and fmax specify the lower and upper freq. for each band, resp.
		fmin = (8., 13.)
		fmax = (13., 30.)
		tmin = .150
		tmax= .350
		sfreq = raw.info['sfreq']  # the sampling frequency

# Now we compute connectivity. To speed things up, we use 2 parallel jobs
# and use mode='fourier', which uses a FFT with a Hanning window
# to compute the spectra (instead of multitaper estimation, which has a
# lower variance but is slower). By using faverage=True, we directly
# average the coherence in the alpha and beta band, i.e., we will only
# get 2 frequency bins
		plv, freqs, times, n_epochs, n_tapers = spectral_connectivity(stcs, method='plv', mode='fourier', indices=indices,sfreq=sfreq, fmin=fmin, fmax=fmax, faverage=True, tmin=tmin, tmax=tmax, n_jobs=2)
#spectral_connectivity(data, method='coh', indices=None, sfreq=6.283185307179586, mode='multitaper', fmin=None, fmax=inf, fskip=0, faverage=False, tmin=None, tmax=None, mt_bandwidth=None, mt_adaptive=False, mt_low_bias=True, cwt_frequencies=None, cwt_n_cycles=7, block_size=1000, n_jobs=1, verbose=None)

		print('Frequencies in Hz over which coherence was averaged for alpha: ')
		print(freqs[0])
		print('Frequencies in Hz over which coherence was averaged for beta: ')
		print(freqs[1])

		# Generate a SourceEstimate with the coherence. This is simple since we
		# used a single seed. For more than one seeds we would have to split coh.
		# Note: We use a hack to save the frequency axis as time
		tmin = np.mean(freqs[0])
		tstep = np.mean(freqs[1]) - tmin
		plv_stc = mne.SourceEstimate(plv, vertices=stc.vertno, tmin=1e-3 * tmin, tstep=1e-3 * tstep, subject=subjects[0])
		data_out = data_path + meg +  'SemLoc_Concrete_PhaseLocking_AlphaBeta_' + label_name[0:-3]
		plv_stc.save(data_out)

		# Now we c+an visualize the coherence using the plot method
		#brain = coh_stc.plot(subjects[0], 'inflated', 'lh', fmin=0.25, fmid=0.4,
		#	             fmax=0.65, time_label='Coherence %0.1f Hz',
		#	             subjects_dir=subjects_dir)
		#brain.show_view('lateral')
