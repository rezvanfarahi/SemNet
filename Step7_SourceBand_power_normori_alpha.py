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
#import matplotlib.pyplot as plt
import numpy as np
import mne
from mne import io
from mne.io import Raw
#import pylab as pl
#from mne import (fiff, read_evokeds, equalize_channels,read_proj, read_selection)
#from mne.preprocessing.ica import ICA


from mne.minimum_norm import (read_inverse_operator, make_inverse_operator, write_inverse_operator, apply_inverse, apply_inverse_epochs, source_induced_power,source_band_induced_power)
from mne.minimum_norm.inverse import (prepare_inverse_operator)

#from surfer import Brain
#from surfer.io import read_stc
#import logging

import scipy.io
#from mne import filter
#from mne.epochs import combine_event_ids
#from mne.layouts import read_layout
from mne.time_frequency import compute_raw_psd
#from mne.connectivity import seed_target_indices, spectral_connectivity

###############################################################################
data_path = '/imaging/rf02/Semnet/' # root directory for your MEG data
#os.chdir(data_path)
subjects_dir = '/imaging/rf02/Semnet/'    # where your MRI subdirectories are
# where event-files are
event_path = '/imaging/rf02/Semnet/'


#label_path = '/group/erp/data/olaf.hauk/MEG/TypLexM/labels/createdlabels_SL/'

inv_path = '/imaging/rf02/Semnet/'
inv_fname = 'InvOp_ico5newreg_fft_pnt1_48_clean_ica_EMEG-inv.fif'#inv_fname = 'InverseOperator_fft_1_48_clean_ica_EMEG-inv.fif'

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
list_all =  ['/meg16_0030/160216/', #0
            '/meg16_0032/160218/', #1
            '/meg16_0034/160219/', #3
            '/meg16_0035/160222/', #4
            '/meg16_0042/160229/', #7
            '/meg16_0045/160303/', #8
            '/meg16_0052/160310/', #10
            '/meg16_0056/160314/',#11
            '/meg16_0069/160405/',#12 
            '/meg16_0070/160407/', #13
            '/meg16_0072/160408/', #15
            '/meg16_0073/160411/', #16
            '/meg16_0075/160411/', #17
            '/meg16_0078/160414/', #18
            '/meg16_0082/160418/', #19
            '/meg16_0086/160422/', #20
            '/meg16_0097/160512/', #21 
            '/meg16_0122/160707/', #22 
            '/meg16_0125/160712/', #24 
            ]



# subjects names used for MRI data
subjects=  ['MRI_meg16_0030' ,#0
            'MRI_meg16_0032' ,#1
            'MRI_meg16_0034' ,#2
            'MRI_meg16_0035' ,#3
            'MRI_meg16_0042' ,#4
            'MRI_meg16_0045' ,#5
            'MRI_meg16_0052' ,#6
            'MRI_meg16_0056' ,#7
    	       'MRI_meg16_0069' ,#8
 	       'MRI_meg16_0070' ,#9
            'MRI_meg16_0072' ,#10
            'MRI_meg16_0073' ,#11
            'MRI_meg16_0075' ,#12
            'MRI_meg16_0078' ,#13
            'MRI_meg16_0082' ,#14
            'MRI_meg16_0086' ,#15
            'MRI_meg16_0097' ,#16
            'MRI_meg16_0122' ,#17
            'MRI_meg16_0125' ,#18
            ]

ll = []
for ss in subject_inds:
 ll.append(list_all[ss])

print "ll:"
print ll



# Labels/ROIs
#labellist = ['atlleft-lh']#, 'atlright-rh'
tmin, tmax = -0.3, 0.6
event_ids = {'visual': 1, 'hear': 2, 'hand': 3, 'neutral': 4, 'emotional': 5,'pwordc': 6}
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
#    src_fname =  data_path + meg + 'forward_5-3L-EMEG-fwd.fif'
    raw_fname = data_path + meg + 'SemDec_blocks_tsss_filt_pnt1_48_ica_raw.fif'
    print raw_fname
 
##events
    event_fname = event_path + meg + 'SemDec_blocks_tsss_filt_pnt1_48_ica_raw-eve.fif'

## print "reading raw data"
    print "reading raw file"
    raw = mne.io.Raw(raw_fname, preload=True)
  
## print "read events"
    print "Reading events from" + event_fname
    events = mne.read_events(event_fname)
    events[:,0] = events[:,0]+np.round( raw.info['sfreq']*stim_delay )

## get configuration info
    include = []
    exclude = raw.info['bads'] # bads

    fname_inv = inv_path + meg + inv_fname
    print "Read inverse operator " + fname_inv
    inv_op = read_inverse_operator(fname_inv)

## epoching
   
    picks = mne.pick_types(raw.info, eeg=True, meg=True, eog=False, exclude='bads')
    print "Epoching"
    if subject_inds[0] in np.array([3,5]):#np.array([4,8,9]):
        reject = dict(eeg=150e-6, grad=200e-12, mag=4e-12)
    else:
        reject = dict(eeg=120e-6, grad=200e-12, mag=4e-12)#, eog=150e-6)
    event_ids = {'visual': 1, 'hear': 2, 'hand': 3, 'pwordc': 6}#'neutral': 4, 'emotional': 5,
    epochs = mne.Epochs(raw, events, event_ids, tmin, tmax, picks=picks, proj=True, baseline=(None, 0), reject=reject,preload=True)
    epochs.resample(200)
#    for eegcnt in range(71):
#        if eegcnt<10:
#            thiseegbad=sum(epochs.drop_log,[]).count('EEG00'+str(eegcnt))
#            if thiseegbad>=90:
#                raw.info['bads'].append(u'EEG00'+str(eegcnt))                
#        else:
#            thiseegbad=sum(epochs.drop_log,[]).count('EEG0'+str(eegcnt))
#            if thiseegbad>=90:
#                raw.info['bads'].append(u'EEG0'+str(eegcnt))
#    picks = mne.pick_types(raw.info, eeg=True, meg=True, eog=False, exclude='bads')
    epochs_list=range(4)
    eid=-1
    for event_id in np.array([1,2,3,6]):#range(6):
        eid=eid+1
#        event_id = eid+1  
        epochs_list[eid]=epochs[epochs.events[:,2]==event_id]
    mne.epochs.equalize_epoch_counts(epochs_list,method='mintime')    
#    epochs_concrete = epochs['cncrt_wrd']
#    epochs_abstract = epochs['abs_wrd']
#    mne.epochs.equalize_epoch_counts([epochs_concrete,epochs_abstract], method='mintime')
    b='alpha'
    
    cwt_frequencies = np.arange(8, 13, 1)#np.array([6,10,22,38])#
    cwt_n_cycles = 4.#cwt_frequencies / float(3)
    print "source_band_induced_power concrete"
    subject_from = subjects[subject_no]
    subject_to = 'fsaverage'
    vertices_to = [np.arange(10242), np.arange(10242)]
    print 'Done epoching'
    event_names=['Visual','Hear','Hand','Pword']
    this_event=-1
    for this_epoch in epochs_list:
        this_event=this_event+1
        epdata=this_epoch.get_data()
        em=np.abs(epdata)
        #emax=np.mean(em[:,:,600:900],axis=2)#
        sa=em[:,:,70:150]#em[:,:,550:950]#.max(axis=2)
        #bm=np.sqrt(np.mean(epdata[:,:,100:500]**2,axis=2))
        bv=np.var(em[:,:,:60],axis=2)#100:500
        edv=sa.copy()
        for ii in range(sa.shape[2]):
            edv[:,:,ii]=(sa[:,:,ii]**2)/bv
        
        ediv=np.mean(np.sqrt(edv),axis=2)
        print ediv.mean()
        snr = ediv.mean()
        lambda2 = 1.0 / snr ** 2
        	#source_induced_power(epochs, inverse_operator, frequencies, label=None, lambda2=0.1111111111111111, method='dSPM', nave=1, n_cycles=5, decim=1, use_fft=False, pick_ori=None, baseline=None, baseline_mode='logratio', pca=True, n_jobs=1, zero_mean=False, verbose=None, pick_normal=None)
        
        	# Compute a source estimate per frequency band
        #bands = dict(theta=[4,7])#, alpha=[8,12], beta=[13, 30], gamma=[31, 45])#, beta=[13, 30], gamma=[31, 45])
        
        concrete_power,concrete_plv=source_induced_power(epochs=this_epoch, inverse_operator=inv_op, frequencies=cwt_frequencies, label=None, lambda2=lambda2, method='MNE', n_cycles=cwt_n_cycles, use_fft=False, baseline=(-0.3,0), baseline_mode='ratio', subject=subject_from,subjects_dir=data_path,n_jobs=4)
        
        #stcs_concrete = source_band_induced_power(epochs_concrete, inv_op, bands=bands, label=None, lambda2=lambda2, method='MNE', n_cycles=2, baseline=(-400,0), baseline_mode='ratio', use_fft=False,n_jobs=4,subject=subject_from,subjects_dir=data_path)
        cpower=concrete_power.mean(axis=1)
        cp_stc = mne.SourceEstimate(cpower, vertices=vertices_to, tmin= -0.3, tstep= 0.005, subject=subject_to)
        
        cplv=concrete_plv.mean(axis=1)
        cplv_stc = mne.SourceEstimate(cplv, vertices=vertices_to, tmin= -0.3, tstep= 0.005, subject=subject_to)
        
        data_out = data_path + meg + 'mineMorphed_ico_SemDec_ica_'+event_names[this_event]+'_Power_normori_ratio_equalized_m500_700_200hz_' + b
        cp_stc.save(data_out)
        data_out = data_path + meg + 'mineMorphed_ico_SemDec_ica_'+event_names[this_event]+'_rPLV_normori_ratio_equalized_m500_700_200hz_' + b
        cplv_stc.save(data_out)
        
        this_epoch_even=this_epoch[::2]
        even_power,even_plv=source_induced_power(epochs=this_epoch_even, inverse_operator=inv_op, frequencies=cwt_frequencies, label=None, lambda2=lambda2, method='MNE', n_cycles=cwt_n_cycles, use_fft=False, baseline=(-0.3,0), baseline_mode='ratio', subject=subject_from,subjects_dir=data_path,n_jobs=4)
        
        #stcs_concrete = source_band_induced_power(epochs_concrete, inv_op, bands=bands, label=None, lambda2=lambda2, method='MNE', n_cycles=2, baseline=(-400,0), baseline_mode='ratio', use_fft=False,n_jobs=4,subject=subject_from,subjects_dir=data_path)
        cpowere=even_power.mean(axis=1)
        cpe_stc = mne.SourceEstimate(cpowere, vertices=vertices_to, tmin= -0.3, tstep= 0.005, subject=subject_to)
        
        cplve=even_plv.mean(axis=1)
        cplve_stc = mne.SourceEstimate(cplve, vertices=vertices_to, tmin= -0.3, tstep= 0.005, subject=subject_to)
        
        data_out = data_path + meg + 'mineMorphed_ico_SemDec_ica_'+event_names[this_event]+'_Power_Even_normori_ratio_equalized_m500_700_200hz_' + b
        cpe_stc.save(data_out)
        data_out = data_path + meg + 'mineMorphed_ico_SemDec_ica_'+event_names[this_event]+'_rPLV_Even_normori_ratio_equalized_m500_700_200hz_' + b
        cplve_stc.save(data_out)
        
        this_epoch_odd=this_epoch[1::2]
        odd_power,odd_plv=source_induced_power(epochs=this_epoch_odd, inverse_operator=inv_op, frequencies=cwt_frequencies, label=None, lambda2=lambda2, method='MNE', n_cycles=cwt_n_cycles, use_fft=False, baseline=(-0.3,0), baseline_mode='ratio', subject=subject_from,subjects_dir=data_path,n_jobs=4)
        
        #stcs_concrete = source_band_induced_power(epochs_concrete, inv_op, bands=bands, label=None, lambda2=lambda2, method='MNE', n_cycles=2, baseline=(-400,0), baseline_mode='ratio', use_fft=False,n_jobs=4,subject=subject_from,subjects_dir=data_path)
        cpowero=odd_power.mean(axis=1)
        cpo_stc = mne.SourceEstimate(cpowero, vertices=vertices_to, tmin= -0.3, tstep= 0.005, subject=subject_to)
        
        cplvo=odd_plv.mean(axis=1)
        cplvo_stc = mne.SourceEstimate(cplvo, vertices=vertices_to, tmin= -0.3, tstep= 0.005, subject=subject_to)
        
        data_out = data_path + meg + 'mineMorphed_ico_SemDec_ica_'+event_names[this_event]+'_Power_Odd_normori_ratio_equalized_m500_700_200hz_' + b
        cpo_stc.save(data_out)
        data_out = data_path + meg + 'mineMorphed_ico_SemDec_ica_'+event_names[this_event]+'_rPLV_Odd_normori_ratio_equalized_m500_700_200hz_' + b
        cplvo_stc.save(data_out)
    
   

#    abstract_power,abstract_plv=source_induced_power(epochs=epochs_abstract, inverse_operator=inv_op, frequencies=cwt_frequencies, label=None, lambda2=lambda2, method='MNE', n_cycles=cwt_n_cycles, use_fft=False, baseline=(-0.4,0), baseline_mode='ratio', subject=subject_from,subjects_dir=data_path,n_jobs=4)
#    
#    #stcs_abstract = source_band_induced_power(epochs_abstract, inv_op, bands=bands, label=None, lambda2=lambda2, method='MNE', n_cycles=2, baseline=(-400,0), baseline_mode='ratio', use_fft=False,n_jobs=4,subject=subject_from,subjects_dir=data_path)
#    apower=abstract_power.mean(axis=1)
#    ap_stc = mne.SourceEstimate(apower, vertices=vertices_to, tmin= -0.5, tstep= 0.005, subject=subject_to)
#    
#    aplv=abstract_plv.mean(axis=1)
#    aplv_stc = mne.SourceEstimate(aplv, vertices=vertices_to, tmin= -0.5, tstep= 0.005, subject=subject_to)
#    
#    data_out = data_path + meg + 'mineMorphed_ico_SemLoc_ica_Abstract_Power2_normori_ratio_equalized_m500_700_200hz_' + b
#    ap_stc.save(data_out)
#    data_out = data_path + meg + 'mineMorphed_ico_SemLoc_ica_Abstract_rPLV2_normori_ratio_equalized_m500_700_200hz_' + b
#    aplv_stc.save(data_out)


    
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

