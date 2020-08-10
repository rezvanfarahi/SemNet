"""
=========================================================
TF representation of TypLex_PW for frequency bands in source labels
=========================================================

"""
# Authors: Alexandre Gramfort <gramfort@nmr.mgh.harvard.edu>
#
# License: BSD (3-clause)


print(__doc__)

import sys
sys.path.append('/imaging/local/software/python_packages/nibabel/1.3.0')
sys.path.append('/imaging/local/software/python_packages/pysurfer/v0.4')
sys.path.append('/imaging/local/software/python_packages/scikit-learn/v0.14.1')
sys.path.append('/imaging/local/software/mne_python/latest')

# for qsub
# (add ! here if needed) 
#sys.path.append('/imaging/local/software/anaconda/latest/x86_64/bin/python')
sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
###
import os
print "hi"
#import matplotlib.pyplot as plt
import numpy as np
print "hi1"
import mne
print "hi2"
from mne import io
print "hi3"
from mne.io import Raw
print "hi4"
#import pylab as pl
#from mne import (fiff, read_evokeds, equalize_channels,read_proj, read_selection)
#from mne.preprocessing.ica import ICA
import scipy.io as sio
print "hi5"
import operator
print "hi6"

from mne.minimum_norm import (read_inverse_operator, make_inverse_operator, write_inverse_operator, apply_inverse, apply_inverse_epochs, source_induced_power,source_band_induced_power)
print "hi7"
from mne.minimum_norm.inverse import (prepare_inverse_operator)
print "hi8"

#from surfer import Brain
#from surfer.io import read_stc
#import logging

import sklearn
print "hi9"
import scipy.io
print "hi10"
#from mne import filter
from mne import find_events
print "hi11"
#from mne.epochs import combine_event_ids
#from mne.layouts import read_layout
from mne.time_frequency import compute_raw_psd, induced_power
print "hi12"
#from mne.connectivity import seed_target_indices, spectral_connectivity

###############################################################################
data_path = '/imaging/rf02/TypLexMEG/' # root directory for your MEG data
#os.chdir(data_path)
subjects_dir = '/imaging/rf02/TypLexMEG/'    # where your MRI subdirectories are
# where event-files are
event_path = '/imaging/rf02/TypLexMEG/'    # where event files are

#label_path = '/group/erp/data/olaf.hauk/MEG/TypLexM/labels/createdlabels_SL/'

inv_path = '/imaging/rf02/TypLexMEG/'
inv_fname = 'TypLex_PW_InverseOperator_EMEG-inv.fif'

# get indices for subjects to be processed from command line input
# 
print sys.argv
subject_inds = []
for ss in sys.argv[1:]:
   subject_inds.append( int( ss ) )

#subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
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
'meg11_0147/110603/'
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



# Labels/ROIs
#labellist = ['atlleft-lh']#, 'atlright-rh'
tmin, tmax = -0.5, 0.7
reject = dict(eeg=120e-6, eog=150e-6, grad=200e-12, mag=4e-12)
# frequency bands with variable number of cycles for wavelets
stim_delay = 0.034 # delay in s
"""
frequencies = np.arange(8, 45, 1)
n_cycles = frequencies / float(7)
n_cycles[frequencies<=20] = 2
n_cycles[frequencies<=20] += np.arange(0,1,0.077)
    # n_cycles[np.logical_and((frequencies>10), (frequencies<=20))] = 3
"""


##loop across subjects...
for ii, meg in enumerate(ll):
    print ii, meg
    if subject_inds[ii] < 3:
    	event_ids = {'TypWords': 2, 'AtypWords': 1, 'TypPWords': 6, 'AtypPWords': 5}
    else:
    	event_ids = {'TypWords': 8, 'AtypWords': 7, 'TypPWords': 6, 'AtypPWords': 5}


## Get raw data
    src_fname =  data_path + meg + 'forward_5-3L-EMEG-fwd.fif'
    raw_fname = data_path + meg + 'typlex_pw_ssstf_fft49_raw.fif'
    print raw_fname
 
##events
    event_fname = event_path + meg + 'typlex_pw_ssstf_fft49_raw-eve.fif'

## print "reading raw data"
    print "reading raw file"
    raw = mne.io.Raw(raw_fname, preload=True)
  
## print "read events"
    print "Reading events from" + event_fname
    events = mne.read_events(event_fname)
    events[:,0] += np.round( raw.info['sfreq']*stim_delay )

## get configuration info
    include = []
    exclude = raw.info['bads'] # bads

    picks = mne.pick_types(raw.info, meg=True, eeg=True, eog=True,
                                stim=False, include=include, exclude=exclude)

    fname_inv = inv_path + meg + inv_fname
    print "Read inverse operator " + fname_inv
    inv_op = read_inverse_operator(fname_inv)

## epoching
    print "mne.Epochs()"
 
    # get epochs from raw
    epochs = mne.Epochs(raw, events, event_ids, tmin, tmax, picks=picks,
                    baseline=(None, 0), reject=reject, preload=True)

    epochs_TypWords = epochs['TypWords']
    epochs_AtypWords = epochs['AtypWords']
    epochs_TypPWords = epochs['TypPWords']
    epochs_AtypPWords = epochs['AtypPWords']

    print 'Done epoching'

    
	#source_induced_power(epochs, inverse_operator, frequencies, label=None, lambda2=0.1111111111111111, method='dSPM', nave=1, n_cycles=5, decim=1, use_fft=False, pick_ori=None, baseline=None, baseline_mode='logratio', pca=True, n_jobs=1, zero_mean=False, verbose=None, pick_normal=None)

	# Compute a source estimate per frequency band
    bands = dict(alpha=[8,12])#,beta=[13, 30], gamma=[31, 45])#, beta=[13, 30], gamma=[31, 45])
    n_cycles=3
    snr = 3.0
    lambda2 = 1.0 / snr ** 2
    subject_from = subjects[subject_inds[0]]
    subject_to = 'average'
    vertices_to = [np.arange(10242), np.arange(10242)]

    print "source_band_induced_power TypWords"
    stcs_TypWords = source_band_induced_power(epochs_TypWords, inv_op, bands=bands, label=None, lambda2=lambda2, method='MNE', n_cycles=n_cycles, baseline=(-400,0), baseline_mode='logratio',use_fft=False, n_jobs=1)

    for b, stc in stcs_TypWords.iteritems():
	stc_to = mne.morph_data(subject_from, subject_to, stc, subjects_dir=data_path, n_jobs=1, grade=vertices_to)
	data_out = data_path + meg + 'Morphed_TypLex_PW_TypWords_Power_logratio_m500_700_' + b
        stc_to.save(data_out)

    print "source_band_induced_power AtypWords"
    stcs_AtypWords = source_band_induced_power(epochs_AtypWords, inv_op, bands=bands, label=None, lambda2=lambda2, method='MNE', n_cycles=n_cycles, baseline=(-400,0), baseline_mode='logratio',use_fft=False, n_jobs=1)
    
    for b, stc in stcs_AtypWords.iteritems():
        stc_to = mne.morph_data(subject_from, subject_to, stc, subjects_dir=data_path, n_jobs=1, grade=vertices_to)
        data_out = data_path + meg + 'Morphed_TypLex_PW_AtypWords_Power_logratio_m500_700_' + b
        stc_to.save(data_out)

    print "source_band_induced_power TypPWords"
    stcs_TypPWords = source_band_induced_power(epochs_TypPWords, inv_op, bands=bands, label=None, lambda2=lambda2, method='MNE', n_cycles=n_cycles, baseline=(-400,0), baseline_mode='logratio',use_fft=False, n_jobs=1)

    for b, stc in stcs_TypPWords.iteritems():
	stc_to = mne.morph_data(subject_from, subject_to, stc, subjects_dir=data_path, n_jobs=1, grade=vertices_to)
	data_out = data_path + meg + 'Morphed_TypLex_PW_TypPWords_Power_logratio_m500_700_' + b
        stc_to.save(data_out)

    print "source_band_induced_power AtypPWords"
    stcs_AtypPWords = source_band_induced_power(epochs_AtypPWords, inv_op, bands=bands, label=None, lambda2=lambda2, method='MNE', n_cycles=n_cycles, baseline=(-400,0), baseline_mode='logratio',use_fft=False, n_jobs=1)
    
    for b, stc in stcs_AtypPWords.iteritems():
        stc_to = mne.morph_data(subject_from, subject_to, stc, subjects_dir=data_path, n_jobs=1, grade=vertices_to)
        data_out = data_path + meg + 'Morphed_TypLex_PW_AtypPWords_Power_logratio_m500_700_' + b
        stc_to.save(data_out)

    
	# Morph using one method (supplying the vertices in fsaverage's source
	# space makes it faster). Note that for any generic subject, you could do:
	#vertices_to = mne.grade_to_vertices(subject_to, grade=5,subjects_dir=data_path)
	# But fsaverage's source space was set up so we can just do this:
	
	

"""
alphafreq=2*np.ones((1,5),dtype=int)
betafreq=4*np.ones((1,18),dtype=int)
aa=np.concatenate((alphafreq,betafreq),axis=1)

	

	###############################################################################
# plot mean power
import matplotlib.pyplot as plt
plt.plot(stcs['alpha'].times, stcs['alpha'].data.mean(axis=0), label='Alpha')
plt.plot(stcs['beta'].times, stcs['beta'].data.mean(axis=0), label='Beta')
plt.xlabel('Time (ms)')
plt.ylabel('Power')
plt.legend()
plt.title('Mean source induced power')
plt.show()
"""

