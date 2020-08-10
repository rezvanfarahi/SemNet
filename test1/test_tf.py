# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 02:59:22 2016

@author: rf02
"""

"""
===================================================
Compute induced power in the source space with dSPM
===================================================

Returns STC files ie source estimates of induced power
for different bands in the source space. The inverse method
is linear based on dSPM inverse operator.

"""
# Authors: Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
#
# License: BSD (3-clause)
import sys
sys.path.append('/imaging/local/software/python_packages/nibabel/1.3.0')
sys.path.append('/imaging/local/software/python_packages/pysurfer/v0.3.1')
sys.path.append('/imaging/local/software/python_packages/scikit-learn/v0.14.1')
sys.path.insert(1,'/imaging/local/software/mne_python/latest_v0.11')

# for qsub
# (add ! here if needed) 
sys.path.append('/imaging/local/software/anaconda/latest/x86_64/bin/python')
sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
###

import numpy as np

import matplotlib.pyplot as plt

import mne
from mne import io
from mne.datasets import sample
from mne.minimum_norm import read_inverse_operator, source_band_induced_power

print(__doc__)

###############################################################################
# Set parameters
data_path = sample.data_path()
raw_fname = data_path + '/MEG/sample/sample_audvis_raw.fif'
fname_inv = data_path + '/MEG/sample/sample_audvis-meg-oct-6-meg-inv.fif'
tmin, tmax, event_id = -0.3, 0.8, 1

# Setup for reading the raw data
raw = io.read_raw_fif(raw_fname)
events = mne.find_events(raw, stim_channel='STI 014')
inverse_operator = read_inverse_operator(fname_inv)

include = []
raw.info['bads'] += ['MEG 2443', 'EEG 053']  # bads + 2 more

# picks MEG gradiometers
picks = mne.pick_types(raw.info, meg=True, eeg=False, eog=True,
                       stim=False, include=include, exclude='bads')

# Load condition 1
event_id = 1
events = events[:10]  # take 10 events to keep the computation time low
# Use linear detrend to reduce any edge artifacts
epochs = mne.Epochs(raw, events, event_id, tmin, tmax, picks=picks,
                    baseline=(None, 0), reject=dict(grad=4000e-13, eog=150e-6),
                    preload=True, detrend=1)

# Compute a source estimate per frequency band
bands = dict(alpha=[8, 12], beta=[12, 30])

stcs = source_band_induced_power(epochs, inverse_operator, bands, n_cycles=5.,
                                 use_fft=False, n_jobs=1)

for b, stc in stcs.iteritems():
    stc.save('induced_power_%s' % b)

###############################################################################
# plot mean power
plt.plot(stcs['alpha'].times, stcs['alpha'].data.mean(axis=0), label='Theta')
plt.plot(stcs['beta'].times, stcs['beta'].data.mean(axis=0), label='Alpha')
plt.xlabel('Time (ms)')
plt.ylabel('Power')
plt.legend()
plt.title('Mean source induced power')
plt.show()