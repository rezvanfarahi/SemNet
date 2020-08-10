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
# (add ! here if needed) 
#sys.path.append('/imaging/local/software/anaconda/latest/x86_64/bin/python')
sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
###
import os
#import matplotlib.pyplot as plt
import numpy as np
import mne
from mne import io
from mne.io import Raw
#import pylab as pl
#from mne import (fiff, read_evokeds, equalize_channels,read_proj, read_selection)
#from mne.preprocessing.ica import ICA
import scipy.io as sio
import operator


#from surfer import Brain
#from surfer.io import read_stc
#import logging

import sklearn
import scipy.io
#from mne import filter
from mne import find_events
#from mne.epochs import combine_event_ids
#from mne.layouts import read_layout

#from mne.connectivity import seed_target_indices, spectral_connectivity

###############################################################################
data_path = '/imaging/rf02/TypLexMEG/' # root directory for your MEG data
orig_path='/group/erp/data/olaf.hauk/MEG/TypLexM/data/MF22_new/'
subjects_dir = '/imaging/rf02/TypLexMEG/'    # where your MRI subdirectories are
# where event-files are
event_path = '/imaging/rf02/TypLexMEG/'    # where event files are

#label_path = '/group/erp/data/olaf.hauk/MEG/TypLexM/labels/createdlabels_SL/'

inv_path = '/imaging/rf02/TypLexMEG/'
inv_fname = 'InverseOperator_fft49_EMEG-inv.fif'

# get indices for subjects to be processed from command line input
# 
print sys.argv
subject_inds = []
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
'meg11_0147/110603/'
]

ll = []
for ss in subject_inds:
 ll.append(list_all[ss])

print "ll:"
print ll


# Labels/ROIs
#labellist = ['atlleft-lh']#, 'atlright-rh'
tmin, tmax = -0.5, 0.7
event_ids = {'cncrt_wrd': 1, 'abs_wrd': 2}
# frequency bands with variable number of cycles for wavelets
stim_delay = 0.034 # delay in s

bb=np.zeros((17,2))
cc=np.zeros((17,2))

##loop across subjects...
for ii,meg in enumerate(ll):
    print ii
## Get raw data
    raw_fname = data_path + meg + 'semloc_ssstf_fft_1_100_clean_ica_raw.fif'
 
##events
    events_fname = orig_path + meg + 'semloc_raw_ssstnew.txt'

## print "reading raw data"
    print "reading raw file"
    raw = mne.io.Raw(raw_fname)

## print "read events"
    print "Reading events"
    events = mne.read_events(events_fname)
    events[:,0] += np.round( raw.info['sfreq']*stim_delay )

## get configuration info
    include = []
    exclude = raw.info['bads'] # bads

    picks = mne.pick_types(raw.info, meg=True, eeg=False, eog=True, stim=False, include=include, exclude=exclude)
## epoching
    print "mne.Epochs()"
 
    # get epochs from raw
    if ii==3:
   	reject = dict( grad=200e-12, mag=4e-12)# eog=150e-6,eeg=150e-6,
    else:
   	reject = dict( grad=200e-12, mag=4e-12)# eog=150e-6,eeg=120e-6,

    epochs = mne.Epochs(raw, events, event_ids, tmin, tmax, picks=picks, baseline=(None, 0), reject=reject, preload=True, detrend=1, add_eeg_ref=True)
    #epochs1 = mne.Epochs(raw1, events, event_ids, tmin, tmax, picks=picks, baseline=(None, 0), reject=reject, preload=True, detrend=1, add_eeg_ref=True)
    
    #epochs.equalize_event_counts(['cncrt_wrd', 'abs_wrd'], copy=False)
    epochs_concrete = epochs['cncrt_wrd']
    epochs_abstract = epochs['abs_wrd']
    aa=epochs_concrete.get_data()
    bb[ii,0]=aa.shape[0]
    aa=epochs_abstract.get_data()
    bb[ii,1]=aa.shape[0]
    print bb[ii,:]
    print 'Done epoching'
    #bb=raw.ch_names[300:400]
    #aa=raw.pick_channels[bb]  
   
# Now let's combine some conditions
"""
# average epochs and get Evoked datasets
    evokeds = [epochs[cond].average() for cond in ['cncrt_wrd', 'abs_wrd']]
    #evokeds1 = [epochs1[cond].average() for cond in ['cncrt_wrd', 'abs_wrd']]
import matplotlib.pyplot as plt
plt.clf()
evokeds[0].plot()
plt.clf()
evokeds1[0].plot()
###############################################################################
# View evoked response

import matplotlib.pyplot as plt
plt.clf()
ax = plt.subplot(2, 1, 1)
evokeds[0].plot(axes=ax)
plt.title('EEG evoked potential, auditory trials')
plt.ylabel('Potential (uV)')
ax = plt.subplot(2, 1, 2)
evokeds[1].plot(axes=ax)
plt.title('EEG evoked potential, visual trials')
plt.ylabel('Potential (uV)')
plt.show()
"""