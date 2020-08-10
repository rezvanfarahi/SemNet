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

#sys.path.append('/imaging/local/software/python_packages/scikit-learn/v0.14.1')
sys.path.append('/imaging/local/software/mne_python/latest_v0.9')

# for qsub
# (add ! here if needed) 
sys.path.append('/imaging/local/software/EPD/latest/x86_64/bin/python')
sys.path.insert(1,'/imaging/local/software/python_packages/scikit-learn/v0.16.1')

#sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
###

import numpy as np
import mne
from mne.io import Raw
from mne.minimum_norm import (read_inverse_operator, apply_inverse, apply_inverse_epochs)
import scipy.io
from mne import find_events
from mne.connectivity import (spectral_connectivity, phase_slope_index)
import numpy as np
from mne import read_labels_from_annot
#import sklearn

data_path = '/imaging/rf02/TypLexMEG/' # root directory for your MEG data
subjects_dir = '/imaging/rf02/TypLexMEG/'    # where your MRI subdirectories are

event_path = event_path = '/group/erp/data/olaf.hauk/MEG/TypLexM/data/MF22_new/'   # where event files are

label_path = '/imaging/rf02/TypLexMEG/fsaverage/label/'

inv_path = '/imaging/rf02/TypLexMEG/'
inv_fname = 'InverseOperator_fft_1_48_clean_ica_EMEG-inv.fif'

# get indices for subjects to be processed from command line input
# 
print sys.argv
subject_inds = []
for ss in sys.argv[1:]:
   subject_inds.append( int( ss ) )



#subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]#0,1,2,3,4,5,6,7,8]#,
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

labels=range(8)
labellist = ['lh.posterior_fusiform','rh.posterior_fusiform', 'lh.lateraloccipital','rh.lateraloccipital', 'lh.ATLsmall','rh.ATLsmall', 'lh.entorhinal', 'rh.entorhinal']#labellist = [ 'leftATL-lh','rightATL-rh']#, 'medialtempright-rh','medialtempleft-lh']'lh.postcentral','rh.postcentral',
for ii,label_name in enumerate(labellist):
    fname_label = label_path + '/' + label_name + '.label'
    labels[ii] = mne.read_label(fname_label)
for meg in ll:        
    print meg
    fname_inv = inv_path + meg + 'InverseOperator_fft_1_48_clean_ica_EMEG-inv.fif'
    fname_raw = data_path + meg + '/semloc_ssstf_fft_1_48_clean_ica_raw.fif'
    fname_event = event_path + meg + 'semloc_raw_ssstnew.txt'

    # Load data
    inverse_operator = read_inverse_operator(fname_inv)
    srcin='/imaging/rf02/TypLexMEG/fsaverage/bem/fsaverage-ico-5-src.fif'
    src = mne.read_source_spaces(srcin)#inverse_operator['src']
    raw = Raw(fname_raw, preload=True)
    
    
    events = mne.read_events(fname_event)
    stim_delay = 0.034
    events[:,0] += np.round( raw.info['sfreq']*stim_delay )
    # Add a bad channel
    #raw.info['bads'] += ['MEG 2443']

    # Pick MEG channels
    picks = mne.pick_types(raw.info, meg=True, eeg=True, stim=False, eog=True, exclude='bads')

    # Define epochs for left-auditory condition
    tmin, tmax = -0.5, 0.7
    if subject_inds[0]==3:
        reject = dict(eeg=150e-6, grad=200e-12, mag=4e-12) #eog=150e-6,
    else:
        reject = dict(eeg=120e-6, grad=200e-12, mag=4e-12) #eog=150e-6, 
    event_ids = {'cncrt_wrd': 1, 'abs_wrd': 2}
    method = "MNE"  # use dSPM method (could also be MNE or sLORETA)

    print "Read epochs"
    epochs = mne.Epochs(raw, events, event_ids, tmin, tmax, picks=picks, baseline=(None, 0), proj=True, reject=reject)
    epochs_concrete = epochs['cncrt_wrd']
    epochs_abstract = epochs['abs_wrd']
    #epochs=epochs_abstract
    mne.epochs.equalize_epoch_counts([epochs_concrete,epochs_abstract],method='mintime')
    events_no_concrete=epochs_concrete.copy().get_data().shape[0]
    events_no_abstract=epochs_abstract.copy().get_data().shape[0]
    stcs_concrete=range(events_no_concrete)
    for ccc in range(len(stcs_concrete)):
        print ccc		
        data_in = data_path + meg + 'Morphed_SemLoc_epochs' + str(ccc) + '_icaclean_Concrete_Source_Evoked_m500_700'
        stcbefore=mne.read_source_estimate(data_in)
        #stcbefore.resample(200)
        stcs_concrete[ccc]=stcbefore
        
    stcs_abstract=range(events_no_abstract)
    for ccc in range(len(stcs_abstract)):
        print ccc		
        data_in = data_path + meg + 'Morphed_SemLoc_epochs' + str(ccc) + '_icaclean_Abstract_Source_Evoked_m500_700'
        stcbefore=mne.read_source_estimate(data_in)
        #stcbefore.resample(200)
        stcs_abstract[ccc]=stcbefore
    method = "MNE"  # use dSPM method (could also be MNE or sLORETA)
    # Get labels for FreeSurfer 'aparc' cortical parcellation with 34 labels/hemi
    label_colors = [label.color for label in labels] 
    label_mne_ctf_mn_ts_concrete = mne.extract_label_time_course(stcs_concrete, labels, src, mode='pca_flip',return_generator=False)
    label_mne_ctf_mn_ts_abstract = mne.extract_label_time_course(stcs_abstract, labels, src, mode='pca_flip',return_generator=False)
    label_mne_ctf_mn_ts_concrete=np.asarray(label_mne_ctf_mn_ts_concrete)
    label_mne_ctf_mn_ts_abstract=np.asarray(label_mne_ctf_mn_ts_abstract)
    
    fmin = (4., 8., 13., 31.)
    fmax = (7., 12., 30., 45.)
    cwt_frequencies = np.arange(4, 46, 1)#np.array([6,10,22,38])#
    cwt_n_cycles = cwt_frequencies / float(2)#np.array([2,3,5,7])#
    sfreq = raw.info['sfreq']  # the sampling frequency
    con_cncrt, freqs1, times1, n_epochs1, n_tapers1 = phase_slope_index(label_mne_ctf_mn_ts_concrete, mode='cwt_morlet', cwt_n_cycles=cwt_n_cycles, cwt_frequencies=cwt_frequencies, sfreq=sfreq, fmin=fmin, fmax=fmax)

    path_c=data_path+ meg+'pca_coh_concrete_multitaper_8regions'
    np.save(path_c,con_cncrt)
    scipy.io.savemat(path_c,{'con_cncrt':con_cncrt})
    print "con finished"
    
#    con_abs=np.zeros((8,8,2,2))
#    for actc,tc in enumerate(range(150,251,100)): #time count
#        tminw=(tc+500)/1000.#tc+500
#        tmaxw=(tc+700)/1000.#tc+700
#        con_abs[:,:,:,actc], freqs1, times1, n_epochs1, n_tapers1 = phase_slope_index(label_mne_ctf_mn_ts_abstract, mode='multitaper', mt_adaptive=True, sfreq=sfreq, fmin=fmin, fmax=fmax,  tmin=tminw, tmax=tmaxw)
    con_abs, freqs1, times1, n_epochs1, n_tapers1 = phase_slope_index(label_mne_ctf_mn_ts_abstract, mode='cwt_morlet', cwt_n_cycles=cwt_n_cycles, cwt_frequencies=cwt_frequencies, sfreq=sfreq, fmin=fmin, fmax=fmax)
    path_a=data_path+ meg+'pca_coh_abstract_multitaper_8regions'
    np.save(path_a,con_abs)
    scipy.io.savemat(path_a,{'con_abs':con_abs})
    print "abs finished"
    
