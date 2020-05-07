"""
=================================================================
Permutation t-test on source data with spatio-temporal clustering
=================================================================

Tests if the evoked response is significantly different between
conditions across subjects (simulated here using one subject's data).
The multiple comparisons problem is addressed with a cluster-level
permutation test across space and time.

"""

# Authors: Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
#          Eric Larson <larson.eric.d@gmail.com>
# License: BSD (3-clause)

print __doc__

import sys
sys.path.insert(1,'/imaging/local/software/mne_python/latest')
sys.path.insert(1,'/imaging/rf02/scikit-learn-0.15.0')
sys.path.insert(1,'/home/rf02/.local/lib/python2.7/site-packages')
#sys.path.insert(1,'/imaging/local/software/anaconda/2.4.1/2/lib/python2.7/site-packages/scipy')
#import scipy
#print scipy.__version__
# for qsub
# (add ! here if needed) /imaging/local/software/anaconda/latest/x86_64/bin/python
#sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
###

import numpy as np
import mne

from mne.minimum_norm import apply_inverse, read_inverse_operator
from mne.epochs import equalize_epoch_counts
from mne import io
import os
from mne.io.pick import channel_indices_by_type


###############################################################################
data_path = '/imaging/rf02/Semnet/' # root directory for your MEG data
orig_path=data_path
subjects_dir = '/imaging/rf02/Semnet/'    # where your MRI subdirectories are
os.chdir('/imaging/rf02/Semnet/')

# get indices for subjects to be processed from command line input
# 
print sys.argv
subject_inds = [0]#
for ss in sys.argv[1:]:
    subject_inds.append( int( ss ) )

#subject_inds=[4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]
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
bands=['theta','alpha','beta','gamma']#'delta',
freqsmin=np.array([4.,8.,12.,30.])#np.array([4,8,30])0.5,
freqsmax=np.array([8.,12.,30.,48.])#np.array([8,12,48])4.,
print "subjects:"
print ll
n_subjects=19
stim_delay = 0.034 # delay in s ((NOTE! not in the example, taken from Olaf's script. Rezvan))
tmin, tmax = -0.3, 0.6
for ii, meg in enumerate(ll):
    print subject_inds[0]
    #cov_path=data_path+meg+'noise_cov'
    #noise_cov=mne.read_cov(cov_path)
    raw_fname = data_path + meg+ 'SemDec_blocks_tsss_filt_pnt1_48_ica_raw.fif' #clean_ssp_
    event_fname = orig_path + meg + 'SemDec_blocks_tsss_filt_pnt1_48_ica_raw-eve.fif'
    subjects_dir = data_path 
    
    tmin = -0.3
    tmax = 0.6  # Use a lower tmax to reduce multiple comparisons
    #   Setup for reading the raw data
    raw = io.Raw(raw_fname, preload=True)
    events = mne.read_events(event_fname)
    for evcnte in range(events.shape[0]-1):
        if events[evcnte,2] in np.array([2]):
            if events[evcnte+1,2] in np.array([79,90]):
                events[evcnte,2]=230
        if events[evcnte,2] in np.array([1]):
            if events[evcnte+1,2] in np.array([71]):
                events[evcnte,2]=230
    stim_delay = 0.034 # delay in s
    events[:,0] = events[:,0]+np.round( raw.info['sfreq']*stim_delay )
    ###############################################################################
    # Read epochs for all channels, removing a bad one
    print "epochs"
    picks = mne.pick_types(raw.info, eeg=True, meg=True, eog=False, exclude='bads')
    print "Epoching"
    if subject_inds[0] in np.array([3,5]):#np.array([4,8,9]):
        reject = dict(eeg=150e-6, grad=200e-12, mag=4e-12)
    else:
        reject = dict(eeg=120e-6, grad=200e-12, mag=4e-12)#, eog=150e-6)
    event_ids = {'visual': 1, 'hear': 2, 'hand': 3, 'pwordc': 6}#'neutral': 4, 'emotional': 5,
#    epochspre = mne.Epochs(raw, events, event_ids, tmin, tmax, picks=picks, proj=True, baseline=(None, 0), reject=None,preload=True)
#    chidx=channel_indices_by_type(epochspre.info)
#    drop_inds=list()
#    edata=epochspre.get_data()
#    for epii in range(edata.shape[0]):
#        for ch_type in ['grad','mag','eeg']:
#            for chjj in chidx[ch_type]:
#                thresh=reject[ch_type]
#                if np.max(edata[epii,chjj,:])-np.min(edata[epii,chjj,:])>thresh:
#                    drop_inds.append(epii)
#    drop_inds=np.unique(np.asarray(drop_inds))
#    print len(epochspre)/900.
    raw_orig=raw.copy()
    for bcnt,b in enumerate(bands):
        print b
        raw=raw_orig.copy()
        freqmin=freqsmin[bcnt].copy();freqmax=freqsmax[bcnt].copy()
#        raw.filter(l_freq=freqmin, h_freq=freqmax,  l_trans_bandwidth=0.1,h_trans_bandwidth=0.1, picks=picks, method='fft',filter_length='auto',phase='zero-double')
#        raw.plot_psd(fmax=50,picks=picks,average=False)
        epochs = mne.Epochs(raw, events, event_ids, tmin, tmax, picks=picks, proj=True, baseline=(None, 0), reject=reject,preload=True)
        epochs.savgol_filter(freqmin)
        epochs.savgol_filter(freqmax)
        epochs.plot_psd(fmax=50.)
#        epochs.drop(drop_inds)
        print len(epochs)/600.#    for eegcnt in range(71):
#        if eegcnt<10:
#            thiseegbad=sum(epochs.drop_log,[]).count('EEG00'+str(eegcnt))
#            if thiseegbad>=90:
#                raw.info['bads'].append(u'EEG00'+str(eegcnt))                
#        else:
#            thiseegbad=sum(epochs.drop_log,[]).count('EEG0'+str(eegcnt))
#            if thiseegbad>=90:
#                raw.info['bads'].append(u'EEG0'+str(eegcnt))
#        epochs_list=range(6)
#        for eid in range(6):
#            event_id = eid+1  
#            epochs_list[eid]=epochs[epochs.events[:,2]==event_id]
#        equalize_epoch_counts(epochs_list,method='mintime')
        
        epochs_list=range(4)
        for eid,event_id in enumerate(list(np.array([1,2,3,6]))):
#                event_id = eid+1  
            epochs_list[eid]=epochs[epochs.events[:,2]==event_id]
        equalize_epoch_counts(epochs_list,method='mintime')

                   
                   
                   
	#    Equalize trial counts to eliminate bias (which would otherwise be
	#    introduced by the abs() performed below)

	###############################################################################
        print "Transform to source space"
    
        fname_inv = data_path + meg + 'InvOp_ico5newreg_fft_pnt1_48_clean_ica_'+b+'_EMEG-inv.fif'
    
        method = "MNE"  # use dSPM method (could also be MNE or sLORETA)
        inverse_operator = read_inverse_operator(fname_inv)
#    sample_vertices = [s['vertno'] for s in inverse_operator['src']]
    #    Let's average and compute inverse, resampling to speed things up
#    evoked=epochs.average()    
#    epdata=evoked.data 
#    em=np.abs(epdata)
#    #emax=np.mean(em[:,:,600:900],axis=2)#
#    sa=em[:,550:950]#.max(axis=2)
#    #bm=np.sqrt(np.mean(epdata[:,:,100:500]**2,axis=2))
#    bv=np.var(em[:,100:500],axis=1)
#    edv=sa.copy()
#    for ii in range(sa.shape[1]):
#        edv[:,ii]=(sa[:,ii]**2)/bv
#
#    ediv=np.sqrt(edv)
#    snr = 3.0#np.mean(ediv)
#    print snr
    
        subject_from = subjects[subject_inds[0]]
        subject_to = 'fsaverage'
        vertices_to = [np.arange(10242), np.arange(10242)]
#    stc_list=range(len(epochs_list))
        event_names = ['Visual', 'Hear', 'Hand', 'Pwordc']#'Neutral', 'Emotional',

#    evoked_list=range(len(epochs_list))
#    condition_list=range(len(epochs_list))
        for evcnt in range(len(epochs_list)):
            this_evoked = epochs_list[evcnt].average()
            
            this_evoked.plot()
#            epdata=this_evoked.data 
#            em=np.abs(epdata)
#            #emax=np.mean(em[:,:,600:900],axis=2)#
#            sa=em[:,350:750]#.max(axis=2)
#            #bm=np.sqrt(np.mean(epdata[:,:,100:500]**2,axis=2))
#            bv=np.var(em[:,50:300],axis=1)
#            edv=sa.copy()
#            for ii in range(sa.shape[1]):
#                edv[:,ii]=(sa[:,ii]**2)/bv    
#            ediv=np.sqrt(edv)
            snr = 3.0#np.mean(ediv)
            lambda2 = 1.0 / snr ** 2
            this_condition = apply_inverse(this_evoked, inverse_operator, lambda2, method,pick_ori="normal")
            this_morphmat=mne.compute_morph_matrix(subject_from, subject_to, this_condition.vertices, vertices_to, subjects_dir=data_path)
            this_stc=this_condition.morph_precomputed(subject_to,vertices_to,this_morphmat, subject_from)
    #        this_stc.resample(200)
            data_out = data_path + meg + 'firstMorphed_epochfilt_SemDec_pnt1_48ica_'+event_names[evcnt]+'_Source_signedEvoked_m300_600_'+b
            this_stc.save(data_out)

    #evoked1 = whiten_evoked(evoked11, noise_cov)
    #evoked1.resample(50)