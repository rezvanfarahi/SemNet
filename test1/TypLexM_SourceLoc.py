"""
=========================================================
TF representation of SemLoc for frequency bands in source labels
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

print "/imaging/local/software/mne_python/git-master-0.8/mne-python/"

import numpy as np
import mne
from mne import fiff, read_proj, read_selection
# from mne.datasets import sample
from mne.fiff import Evoked
from mne.minimum_norm import apply_inverse_epochs, apply_inverse, read_inverse_operator, source_induced_power
import pylab as pl
from mne.layouts import read_layout
from mne.time_frequency import compute_raw_psd, induced_power
from mne.connectivity import seed_target_indices, spectral_connectivity

###############################################################################
data_path = '/home/rf02/rezvan/TypLexMEG'	# where subdirs for MEG data are

# where event-files are
event_path = '/home/rf02/rezvan/TypLexMEG'     # where event files are

label_path = '/group/erp/data/olaf.hauk/MEG/TypLexM/labels/createdlabels_SL/'

inv_path = '/group/erp/data/olaf.hauk/MEG/TypLexM/data/MF22_new/'
inv_fname = 'InvOp_3L-EMEG-loose0.2-inv.fif'

# get indices for subjects to be processed from command line input
# 
print sys.argv
subject_inds = []
for ss in sys.argv[1:]:
 subject_inds.append( int( ss ) )


print "subject_inds:"
print subject_inds


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

ll = []
for ss in subject_inds:
 ll.append(list_all[ss])

print "ll:"
print ll


# Labels/ROIs
labellist = ['atlleft-lh', 'atlright-rh']

stim_delay = 0.034 # delay in s

##loop across subjects...
for meg in ll:

## Get raw data
    raw_fname = data_path + meg + 'semloc_raw_ssstf_raw.fif'
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

    picks = fiff.pick_types(raw.info, meg=True, eeg=True, eog=True,
                                stim=False, include=include, exclude=exclude)

    fname_inv = inv_path + meg + inv_fname
    print "Read inverse operator " + fname_inv
    inv_op = read_inverse_operator(fname_inv)

## epoching
    print "mne.Epochs()"
 
    tmin, tmax = -0.2, 0.5
    reject = dict(eeg=120e-6, eog=150e-6, grad=200e-12, mag=4e-12)
    event_ids = {'C': 1, 'A': 2}

    # frequency bands with variable number of cycles for wavelets
    frequencies = np.arange(8, 45, 1)
    n_cycles = frequencies / float(7)
    n_cycles[frequencies<=20] = 2
    n_cycles[frequencies<=20] += np.arange(0,1,0.077)
    # n_cycles[np.logical_and((frequencies>10), (frequencies<=20))] = 3

    # get epochs from raw
    epochs = mne.Epochs(raw, events, event_ids, tmin, tmax, picks=picks,
                    baseline=(None, 0), reject=reject, preload=True)

    epochs_A = epochs['A']
    epochs_C = epochs['C']

    print 'Done epoching'

    for label in labellist:

        print label

        fname_label = label_path + subjects[list_all.index(meg)] + '/' + label + '.label'
        labeldat = mne.read_label(fname_label)
        
        # Source Induced Power and phase locking in label
        power_A, phase_lock_A = source_induced_power(epochs_A, inv_op, frequencies, labeldat, baseline=(-0.2,0), baseline_mode='percent', n_cycles=n_cycles, pca=True, n_jobs=1, method='MNE')
        power_A = np.mean(power_A, axis=0)
        phase_lock_A = np.mean(phase_lock_A, axis=0)

        power_C, phase_lock_C = source_induced_power(epochs_C, inv_op, frequencies, labeldat, baseline=(-0.2,0), baseline_mode='percent', n_cycles=n_cycles, pca=True, n_jobs=1, method='MNE')
        power_C = np.mean(power_C, axis=0)
        phase_lock_C = np.mean(phase_lock_C, axis=0)

        # write source band results for individual subjects
        out_fname = data_path + meg + '/semloc_power_A_' + label + '.np'
        print out_fname
        np.save(out_fname, power_A)
        out_fname = data_path + meg + '/semloc_power_C_' + label + '.np'
        print out_fname
        np.save(out_fname, power_C)
        out_fname = data_path + meg + '/semloc_phase_lock_A_' + label + '.np'
        print out_fname
        np.save(out_fname, phase_lock_A)
        out_fname = data_path + meg + '/semloc_phase_lock_C_' + label + '.np'
        print out_fname
        np.save(out_fname, phase_lock_C)
