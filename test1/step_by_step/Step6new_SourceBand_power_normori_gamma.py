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
sys.path.append('/imaging/local/software/python_packages/scikit-learn/v0.14.1')
#sys.path.insert(1,'/imaging/local/software/mne_python/latest_v0.9')#/imaging/rf02/TypLexMEG/mne_python0.9
sys.path.insert(1,'/imaging/rf02/TypLexMEG/mne_python0.9')
sys.path.insert(1,'/imaging/local/software/python')
sys.path.insert(1,'/imaging/local/software/python_packages/scikit-learn/v0.13.1/lib.linux-x86_64-2.7/sklearn/externals')
import joblib

# for qsub
# (add ! here if needed) 
#sys.path.append('/imaging/local/software/anaconda/latest/x86_64/bin/python')
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


from mne.minimum_norm import (read_inverse_operator, make_inverse_operator, write_inverse_operator, apply_inverse, apply_inverse_epochs, source_induced_power,source_band_induced_power)
print "hi7"
from mne.minimum_norm.inverse import (prepare_inverse_operator)
print "hi8"

#from surfer import Brain
#from surfer.io import read_stc
#import logging

print "hi9"
import scipy.io
#from mne import filter
#from mne.epochs import combine_event_ids
#from mne.layouts import read_layout
from mne.time_frequency import compute_raw_psd
print "hi12"
#from mne.connectivity import seed_target_indices, spectral_connectivity

###############################################################################
data_path = '/imaging/rf02/TypLexMEG/' # root directory for your MEG data
#os.chdir(data_path)
subjects_dir = '/imaging/rf02/TypLexMEG/'    # where your MRI subdirectories are
# where event-files are
event_path = '/group/erp/data/olaf.hauk/MEG/TypLexM/data/MF22_new/'


#label_path = '/group/erp/data/olaf.hauk/MEG/TypLexM/labels/createdlabels_SL/'

inv_path = '/imaging/rf02/TypLexMEG/'
inv_fname = 'InvOp_ico5diagreg_fft_1_48_clean_ica_EMEG-inv.fif'#inv_fname = 'InverseOperator_fft_1_48_clean_ica_EMEG-inv.fif'

# get indices for subjects to be processed from command line input
# 
print sys.argv
subject_inds = []
for ss in sys.argv[1:]:
   subject_inds.append( int( ss ) )

#subject_inds=[6,7,8,9,10,11,12,13,14,15,16]
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



# Labels/ROIs
#labellist = ['atlleft-lh']#, 'atlright-rh'
tmin, tmax = -0.5, 0.7
event_ids = {'cncrt_wrd': 1, 'abs_wrd': 2}
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
for meg in ll:
    subject_no=subject_inds[0]
## Get raw data
    src_fname =  data_path + meg + 'forward_5-3L-EMEG-fwd.fif'
    raw_fname = data_path + meg + 'semloc_ssstf_fft_1_48_clean_ica_raw.fif'
    print raw_fname
 
##events
    event_fname = event_path + meg + 'semloc_raw_ssstnew.txt'

## print "reading raw data"
    print "reading raw file"
    raw = mne.io.Raw(raw_fname, preload=False)
  
## print "read events"
    print "Reading events from" + event_fname
    events = mne.read_events(event_fname)
    events[:,0] = events[:,0]+np.round( raw.info['sfreq']*stim_delay )

## get configuration info
    include = []
    exclude = raw.info['bads'] # bads

    picks = mne.pick_types(raw.info, meg=True, eeg=True, eog=True, stim=False, include=include, exclude=exclude)

    fname_inv = inv_path + meg + inv_fname
    print "Read inverse operator " + fname_inv
    inv_op = read_inverse_operator(fname_inv)

## epoching
    print "mne.Epochs()"
    if subject_no==3:
    	reject = dict(eeg=150e-6, grad=200e-12, mag=4e-12) #eog=150e-6,
    else:
    	reject = dict(eeg=120e-6, grad=200e-12, mag=4e-12) #eog=150e-6,

    # get epochs from raw
    epochs = mne.Epochs(raw, events, event_ids, tmin, tmax, picks=picks, baseline=(None, 0), reject=reject, proj=True, preload=True)
    epochs.resample(200)
    epochs_concrete = epochs['cncrt_wrd']
    epochs_abstract = epochs['abs_wrd']
    mne.epochs.equalize_epoch_counts([epochs_concrete,epochs_abstract], method='mintime')

    print 'Done epoching'

    epdata=epochs.get_data()
    em=np.abs(epdata)
    #emax=np.mean(em[:,:,600:900],axis=2)#
    sa=em[:,:,110:190]#em[:,:,550:950]#.max(axis=2)
    #bm=np.sqrt(np.mean(epdata[:,:,100:500]**2,axis=2))
    bv=np.var(em[:,:,20:100],axis=2)#100:500
    edv=sa.copy()
    for ii in range(sa.shape[2]):
        edv[:,:,ii]=(sa[:,:,ii]**2)/bv
    
    ediv=np.mean(np.sqrt(edv),axis=2)
    print ediv.mean()
    
	#source_induced_power(epochs, inverse_operator, frequencies, label=None, lambda2=0.1111111111111111, method='dSPM', nave=1, n_cycles=5, decim=1, use_fft=False, pick_ori=None, baseline=None, baseline_mode='logratio', pca=True, n_jobs=1, zero_mean=False, verbose=None, pick_normal=None)

	# Compute a source estimate per frequency band
    #bands = dict(theta=[4,7])#, alpha=[8,12], beta=[13, 30], gamma=[31, 45])#, beta=[13, 30], gamma=[31, 45])
    b='gamma'
    snr = ediv.mean()
    lambda2 = 1.0 / snr ** 2
    cwt_frequencies = np.arange(31, 46, 1)#np.array([6,10,22,38])#
    cwt_n_cycles = 9.#cwt_frequencies / float(3)
    print "source_band_induced_power concrete"
    subject_from = subjects[subject_no]
    subject_to = 'fsaverage'
    vertices_to = [np.arange(10242), np.arange(10242)]
    concrete_power,concrete_plv=source_induced_power(epochs=epochs_concrete, inverse_operator=inv_op, frequencies=cwt_frequencies, label=None, lambda2=lambda2, method='MNE', n_cycles=cwt_n_cycles, use_fft=False, baseline=(-0.4,0), baseline_mode='ratio', subject=subject_from,subjects_dir=data_path,n_jobs=4)
    
    #stcs_concrete = source_band_induced_power(epochs_concrete, inv_op, bands=bands, label=None, lambda2=lambda2, method='MNE', n_cycles=2, baseline=(-400,0), baseline_mode='ratio', use_fft=False,n_jobs=4,subject=subject_from,subjects_dir=data_path)
    cpower=concrete_power.mean(axis=1)
    cp_stc = mne.SourceEstimate(cpower, vertices=vertices_to, tmin= -0.5, tstep= 0.005, subject=subject_to)
    
    cplv=concrete_plv.mean(axis=1)
    cplv_stc = mne.SourceEstimate(cplv, vertices=vertices_to, tmin= -0.5, tstep= 0.005, subject=subject_to)
    
    data_out = data_path + meg + 'mineMorphed_ico_SemLoc_ica_Concrete_Power2_normori_ratio_equalized_m500_700_200hz_' + b
    cp_stc.save(data_out)
    data_out = data_path + meg + 'mineMorphed_ico_SemLoc_ica_Concrete_rPLV2_normori_ratio_equalized_m500_700_200hz_' + b
    cplv_stc.save(data_out)
    
   

    abstract_power,abstract_plv=source_induced_power(epochs=epochs_abstract, inverse_operator=inv_op, frequencies=cwt_frequencies, label=None, lambda2=lambda2, method='MNE', n_cycles=cwt_n_cycles, use_fft=False, baseline=(-0.4,0), baseline_mode='ratio', subject=subject_from,subjects_dir=data_path,n_jobs=4)
    
    #stcs_abstract = source_band_induced_power(epochs_abstract, inv_op, bands=bands, label=None, lambda2=lambda2, method='MNE', n_cycles=2, baseline=(-400,0), baseline_mode='ratio', use_fft=False,n_jobs=4,subject=subject_from,subjects_dir=data_path)
    apower=abstract_power.mean(axis=1)
    ap_stc = mne.SourceEstimate(apower, vertices=vertices_to, tmin= -0.5, tstep= 0.005, subject=subject_to)
    
    aplv=abstract_plv.mean(axis=1)
    aplv_stc = mne.SourceEstimate(aplv, vertices=vertices_to, tmin= -0.5, tstep= 0.005, subject=subject_to)
    
    data_out = data_path + meg + 'mineMorphed_ico_SemLoc_ica_Abstract_Power2_normori_ratio_equalized_m500_700_200hz_' + b
    ap_stc.save(data_out)
    data_out = data_path + meg + 'mineMorphed_ico_SemLoc_ica_Abstract_rPLV2_normori_ratio_equalized_m500_700_200hz_' + b
    aplv_stc.save(data_out)


    
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

