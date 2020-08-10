"""
=========================================================
TF representation of SemLoc for frequency bands in source labels
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
# (add ! here if needed) /imaging/local/software/anaconda/latest/x86_64/bin/python
sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
###

import matplotlib.pyplot as plt
import numpy as np
import mne
from mne import io
from mne.io import Raw
import pylab as pl
from mne import (fiff, read_evokeds, equalize_channels,read_proj, read_selection)
from mne.preprocessing.ica import ICA
import scipy.io as sio
import operator

from mne.minimum_norm import (read_inverse_operator, make_inverse_operator, write_inverse_operator, apply_inverse, read_inverse_operator,apply_inverse_epochs, apply_inverse, source_induced_power,source_band_induced_power)
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

from mne.fiff import Evoked
from mne.layouts import read_layout
from mne.time_frequency import compute_raw_psd, induced_power
from mne.connectivity import seed_target_indices, spectral_connectivity

###############################################################################
data_path = '/home/rf02/rezvan/TypLexMEG/' # root directory for your MEG data
subjects_dir = '/home/rf02/rezvan/TypLexMEG/'    # where your MRI subdirectories are
# where event-files are
event_path = '/home/rf02/rezvan/TypLexMEG/'    # where event files are

label_path = '/group/erp/data/olaf.hauk/MEG/TypLexM/labels/createdlabels_SL/'

inv_path = '/home/rf02/rezvan/TypLexMEG/'
inv_fname = 'InverseOperator_EMEG-inv.fif'

# get indices for subjects to be processed from command line input
# 
#print sys.argv
#subject_inds = []
#for ss in sys.argv[1:]:
 #  subject_inds.append( int( ss ) )


#print "subject_inds:"
#print subject_inds


list_all = ['meg10_0378/101209/', 
#'meg11_0050/110307/', 
#'meg11_0052/110307/', 
#'meg11_0069/110315/', 
#'meg11_0086/110322/', 
#'meg11_0091/110328/', 
#'meg11_0096/110404/', 
#'meg11_0101/110411/', 
#'meg11_0102/110411/', 
#'meg11_0104/110412/', 
#'meg11_0112/110505/', 
#'meg11_0118/110509/', 
#'meg11_0131/110519/', 
#'meg11_0144/110602/', 
#'meg11_0147/110603/', 
#'meg11_0026/110223/', 
#'meg10_0390/101214/'
]

# subjects names used for MRI data
subjects=['SS1_Meg10_0378',
#'SS2_Meg10_0390',
#'SS3_Meg11_0026',
#'SS4_Meg11_0050',
#'SS5_Meg11_0052',
#'SS6_Meg11_0069',
#'SS9_Meg11_0086',
#'SS10_Meg11_0091',
#'SS12_Meg11_0096',
#'SS13_Meg11_0101',
#'SS14_Meg11_0102',
#'SS15_Meg11_0112',
#'SS16_Meg11_0104',
#'SS18_Meg11_0118',
#'SS19_Meg11_0131',
#'SS20_Meg11_0144',
#'SS21_Meg11_0147'
]

#ll = []
#for ss in subject_inds:
#   ll.append(list_all[ss])

#print "ll:"
#print ll


# Labels/ROIs
labellist = ['atlleft-lh', 'atlright-rh']
tmin, tmax = -0.2, 0.5
reject = dict(eeg=120e-6, eog=150e-6, grad=200e-12, mag=4e-12)
event_ids = {'abs_wrd': 1, 'cncrt_wrd': 2}
# frequency bands with variable number of cycles for wavelets
stim_delay = 0.034 # delay in s
frequencies = np.arange(8, 45, 1)
n_cycles = frequencies / float(7)
n_cycles[frequencies<=20] = 2
n_cycles[frequencies<=20] += np.arange(0,1,0.077)
    # n_cycles[np.logical_and((frequencies>10), (frequencies<=20))] = 3



##loop across subjects...
for meg in list_all:

## Get raw data
    raw_fname = data_path + meg + '/semloc_raw_ssstf_raw.fif'
    print raw_fname
 
##events
    event_fname = event_path + meg + 'semloc_raw_ssstf-eve.fif'

## print "reading raw data"
    print "reading raw file"
    raw = mne.io.Raw(raw_fname, preload=False)
  
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

    epochs_concrete = epochs['cncrt_wrd']
    epochs_abstract = epochs['abs_wrd']

    print 'Done epoching'

    for label in labellist:

        print label

        fname_label = label_path + subjects[list_all.index(meg)] + '/' + label + '.label'
        labeldat = mne.read_label(fname_label)

       
# subtract the evoked response in order to exclude evoked activity
	epochs_induced_concrete = epochs_concrete.copy().subtract_evoked()
	epochs_induced_abstract = epochs_abstract.copy().subtract_evoked()

	import matplotlib.pyplot as plt
	plt.close('all')

	for ii, (this_epochs, title) in enumerate(zip([epochs_concrete, epochs_induced_concrete], ['evoked + induced', 'induced only'])):
	    print ii
	  # Source Induced Power and phase locking in label
            #power_concrete, phase_lock_concrete = source_induced_power(this_epochs, inv_op, frequencies, label=labeldat, baseline=(-0.2,0.), baseline_mode='percent', n_cycles=n_cycles, pca=True, n_jobs=1, method='MNE')
		
	#source_band_induced_power(epochs, inverse_operator, bands, label=None, lambda2=0.1111111111111111, method='dSPM', nave=1, n_cycles=5, df=1, use_fft=False, decim=1, baseline=None, baseline_mode='logratio', pca=True, n_jobs=1, verbose=None)  
	#source_induced_power(epochs, inverse_operator, frequencies, label=None, lambda2=0.1111111111111111, method='dSPM', nave=1, n_cycles=5, decim=1, use_fft=False, pick_ori=None, baseline=None, baseline_mode='logratio', pca=True, n_jobs=1, zero_mean=False, verbose=None, pick_normal=None)

	# Compute a source estimate per frequency band
	    bands = dict(alpha=[8, 12], beta=[13, 30])
	    snr = 3.0
	    lambda2 = 1.0 / snr ** 2
    	    stcs = source_band_induced_power(this_epochs, inv_op, bands, label=labeldat, lambda2=lambda2, method='MNE', n_cycles=1, use_fft=False, n_jobs=1)

	    for b, stc in stcs.iteritems():
	       stc.save('induced_power_%s_%s_%d' % (label, b, ii))

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

	#for ii, (this_epochs2, title) in enumerate(zip([epochs_abstract, epochs_induced_abstract], ['evoked + induced', 'induced only'])):
	  # Source Induced Power and phase locking in label
         #   power_abstract, phase_lock_abstract = source_induced_power(this_epochs2, inv_op, frequencies, labeldat, baseline=(-0.2,0.), baseline_mode='percent', n_cycles=n_cycles, pca=True, n_jobs=1, method='MNE')

