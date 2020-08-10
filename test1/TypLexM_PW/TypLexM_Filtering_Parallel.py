"""
=========================================================
Test script to compute  evoked data
=========================================================

"""
# Authors: Alexandre Gramfort <gramfort@nmr.mgh.harvard.edu>
#
# License: BSD (3-clause)

print __doc__

# Russell's addition
import sys
sys.path.append('/imaging/local/software/python_packages/nibabel/v1.3.0')
sys.path.append('/imaging/local/software/python_packages/pysurfer/v0.3.1')
# End

# for qsub
# (add ! here if needed) /imaging/local/software/anaconda/latest/x86_64/bin/python
sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
###

import os
import numpy as np
import mne
from mne import io, read_proj, read_selection, filter
# from mne.datasets import sample
from mne.fiff import Evoked
from mne.io import Raw
from mne.minimum_norm import apply_inverse_epochs, apply_inverse, read_inverse_operator, source_induced_power
import pylab as pl
from mne.layouts import read_layout
from mne.time_frequency import compute_raw_psd, induced_power
from mne.connectivity import seed_target_indices, spectral_connectivity
import matplotlib.pyplot as plt
#### mne_check_eeg_locations --file typlex_pw_raw_ssst.fif
###############################################################################
main_path = '/group/erp/data/olaf.hauk/MEG/TypLexM/data/MF22_new/'
data_path = '/imaging/rf02/TypLexMEG/'	# where subdirs for MEG data are
#os.chdir('/home/rf02/rezvan/test1')

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
'meg11_0050/110307/', 
'meg11_0052/110307/', 
#'meg11_0069/110315/', 
'meg11_0086/110322/', 
'meg11_0091/110328/', 
'meg11_0096/110404/', 
'meg11_0101/110411/', 
'meg11_0102/110411/', 
'meg11_0104/110412/', 
'meg11_0112/110505/', 
'meg11_0118/110509/', 
'meg11_0131/110519/', 
'meg11_0144/110602/', 
'meg11_0147/110603/', 
'meg11_0026/110223/', 
]


ll = []
for ss in subject_inds:
 ll.append(list_all[ss])

print "ll:"
print ll

stim_delay = 0.034 # delay in s

##loop across subjects...
for ii, meg in enumerate(ll):
	print ii
	print meg

## Get raw data
	raw_fname = main_path + meg + 'typlex_pw_raw_ssst.fif'#'typlex_pw_raw_ssstf_raw.fif'#
	raw = mne.io.Raw(raw_fname, preload=True)
	include = []
	exclude = raw.info['bads'] # bads
	#r=raw.copy()
	
	picks = mne.pick_types(raw.info, meg=True, eeg=True, eog=True, stim=False, include=include, exclude=exclude)
	raw.filter(l_freq=0.5, h_freq=49, l_trans_bandwidth=0.1, picks=picks, method='fft')
	
	print raw_fname
	
	#raw.plot_psds(area_mode='range', tmax=100, fmin=0, fmax=100, n_fft=2048, n_jobs=1)
	#r.plot_psds(area_mode='range', tmax=100, fmin=0, fmax=100, n_fft=2048, n_jobs=1)
		 
##events
	print "finding filtering"
	#raw=mne.filter.band_pass_filter(raw, Fs=raw.info['sfreq'], Fp1=1, Fp2=49, picks=picks, filter_length='10s', l_trans_bandwidth=0.4, h_trans_bandwidth=0.4, n_jobs=1, method='fft', iir_params=None, verbose=None)

# Writing events
	print "writing filtred"
	out_name=data_path + meg + 'typlex_pw_ssstf_fft49_raw.fif'
#cd out_path
	raw.save(out_name, overwrite=True)
	
