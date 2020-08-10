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
from mne import viz, read_proj, read_selection, io
# from mne.datasets import sample
from mne.fiff import Evoked
from mne.io import Raw
from mne.minimum_norm import apply_inverse_epochs, apply_inverse, read_inverse_operator, source_induced_power
import pylab as pl
from mne.layouts import read_layout
from mne.time_frequency import compute_raw_psd, induced_power
from mne.connectivity import seed_target_indices, spectral_connectivity

###############################################################################
data_path = '/imaging/rf02/TypLexMEG/'	# where subdirs for MEG data are
os.chdir('/imaging/rf02/TypLexMEG/')

# get indices for subjects to be processed from command line input
# 
print sys.argv
subject_inds = [0]
for ss in sys.argv[1:]:
 subject_inds.append( int( ss ) )

#subject_inds=[0,1]
print "subject_inds:"
print subject_inds
print "No rejection"

list_all = ['meg10_0378/101209/', 
'meg10_0390/101214/',
'meg11_0050/110307/', 
'meg11_0052/110307/', 
'meg11_0069/110315/', 
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
for meg in ll:

## Get raw data
	raw_fname = data_path + meg + 'semloc_ssstf_fft49_raw.fif'
	raw = Raw(raw_fname)
	print raw_fname

	mne.viz.plot_raw_psds(raw, tmin=0.0, tmax=800.0, fmin=0, fmax=60, proj=False, n_fft=2048, picks=None, ax=None, color='black', area_mode='std', area_alpha=0.33, n_jobs=1, verbose=None)
