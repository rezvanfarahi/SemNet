# Author: Martin Luessi <mluessi@nmr.mgh.harvard.edu>
#
# License: BSD (3-clause)

print(__doc__)
print "hi! this script is for spectral functional connectivity. Very nice maps are fruits which worth waiting for I bet!"
import sys

#sys.path.append('/imaging/local/software/python_packages/nibabel/1.3.0')
#sys.path.append('/imaging/local/software/python_packages/pysurfer/v0.4')
#sys.path.append('/imaging/local/software/python_packages/scikit-learn/v0.14.1')
#sys.path.append('/imaging/local/software/mne_python/latest')

# for qsub
# (add ! here if needed)
sys.path.append('/imaging/local/software/EPD/latest/x86_64/bin/python')
sys.path.append('/imaging/local/software/mne_python/latest_v0.9')
###


import os
import numpy as np
import mne
from mne.io import Raw
from mne.minimum_norm import (read_inverse_operator, apply_inverse, apply_inverse_epochs)
from mne.connectivity import seed_target_indices, spectral_connectivity
data_path = '/imaging/rf02/TypLexMEG/' # root directory for your MEG data
os.chdir(data_path)
event_path = '/group/erp/data/olaf.hauk/MEG/TypLexM/data/MF22_new/'   # where event files are
label_path = '/imaging/rf02/TypLexMEG/fsaverage/label'
inv_path = '/imaging/rf02/TypLexMEG/'
subjects_dir = data_path
print sys.argv
subject_inds = []
for ss in sys.argv[1:]:
   subject_inds.append( int( ss ) )

#subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
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




#label_name_lh = ['atlleft-lh', 'atlright-rh']


# Labels/ROIs
labellist = [ 'lh.precentral','rh.precentral','lh.postcentral','rh.postcentral', 'lh.lateraloccipital','rh.lateraloccipital']#labellist = [ 'posterior_fusiform-lh','rh.lateralorbitofrontal', 'rh.medialorbitofrontal']#labellist = [ 'leftATL-lh','rightATL-rh']#, 'medialtempright-rh','medialtempleft-lh']

ii=-1
for meg in ll:
    print subject_inds
    tmin, tmax = -0.5, 0.7
    if subject_inds[0]==3:
        reject = dict(eeg=150e-6, grad=200e-12, mag=4e-12) #eog=150e-6,
    else:
        reject = dict(eeg=120e-6, grad=200e-12, mag=4e-12) #eog=150e-6,
    event_ids = {'cncrt_wrd': 1, 'abs_wrd': 2}
# frequency bands with variable number of cycles for wavelets
    stim_delay = 0.034 # delay in s
#  n_cycles[np.logical_and((frequencies>10), (frequencies<=20))] = 3

    method = "MNE"  # use dSPM method (could also be MNE or sLORETA)
    ii=ii+1
    fname_inv = inv_path + meg + 'InverseOperator_fft_1_48_clean_ica_EMEG-inv.fif'
    fname_raw = data_path + meg + '/semloc_ssstf_fft_1_48_clean_ica_raw.fif'
    fname_event = event_path + meg + 'semloc_raw_ssstnew.txt'
# Load data
    inverse_operator = read_inverse_operator(fname_inv)
    raw = Raw(fname_raw)
    events = mne.read_events(fname_event)
# pick MEG channels
    picks = mne.pick_types(raw.info, meg=True, eeg=True, stim=False, eog=True, exclude='bads')

# Read epochs
    epochs = mne.Epochs(raw, events, event_ids, tmin, tmax, picks=picks, baseline=(None, 0), proj=True, reject=reject)
    epochs_concrete = epochs['cncrt_wrd']
    epochs_abstract = epochs['abs_wrd']

    mne.epochs.equalize_epoch_counts([epochs_concrete,epochs_abstract])
    epochs=epochs_concrete
    events_no=epochs.copy().get_data().shape[0]
    stcs=range(events_no)
    for ccc in range(len(stcs)):
        print ccc		
        data_in = data_path + meg + 'Morphed_SemLoc_epochs' + str(ccc) + '_icaclean_Concrete_Source_Evoked_m500_700'
        stcs[ccc]=mne.read_source_estimate(data_in)

# First, we find the most active vertex in the left auditory cortex, which
#  we will later use as seed for the connectivity computation
    snr = 3.0
    lambda2 = 1.0 / snr ** 2
    evoked = epochs.average()
    sct_incoh = data_path + meg + 'firstMorphed_SemLoc_icaclean_Concrete_Source_Evoked_m500_700'
    stc_coh = mne.read_source_estimate(sct_incoh)

    for label_name in labellist:
        print label_name
        tmin, tmax = -0.5, 0.7
        if subject_inds[0]==3:
            reject = dict(eeg=150e-6, grad=200e-12, mag=4e-12) #eog=150e-6,
        else:
            reject = dict(eeg=120e-6, grad=200e-12, mag=4e-12) #eog=150e-6,
        event_ids = {'cncrt_wrd': 1, 'abs_wrd': 2}
# frequency bands with variable number of cycles for wavelets
        frequencies = np.arange(8, 45, 1)
        n_cycles = frequencies / float(7)
        n_cycles[frequencies<=20] = 2
        n_cycles[frequencies<=20] += np.arange(0,1,0.077)
# n_cycles[np.logical_and((frequencies>10), (frequencies<=20))] = 3
        subject_from = subjects[subject_inds[0]]
        subject_to = 'fsaverage'
        vertices_to = [np.arange(10242), np.arange(10242)]

        method = "MNE"  # use dSPM method (could also be MNE or sLORETA)

        fname_label = label_path + '/' + label_name + '.label'
        label1 = mne.read_label(fname_label)
# Restrict the source estimate to the label in the left auditory cortex
        #stcp = apply_inverse(evoked, inverse_operator, lambda2, method, pick_ori="normal")
        #stc = mne.morph_data(subject_from, subject_to, stcp, subjects_dir=data_path, n_jobs=1, grade=vertices_to)
        stc_label = stc_coh.in_label(label1)

# Find number and index of vertex with most power
        src_pow = np.sum(stc_label.data ** 2, axis=1)
        if label_name[0] is 'l':
            print label_name
            seed_vertno = stc_label.vertices[0][np.argmax(src_pow)]
            seed_idx = np.searchsorted(stc_coh.vertices[0], seed_vertno)  # index in original stc
            print seed_idx

        else:
            print label_name
            seed_vertno = stc_label.vertices[1][np.argmax(src_pow)]
            seed_idx = stc_coh.vertices[0].shape[0]+ np.searchsorted(stc_coh.vertices[1], seed_vertno)  # index in original stc
            print seed_idx

# Generate index parameter for seed-based connectivity analysis
        n_sources = stc_coh.data.shape[0]
        indices = seed_target_indices([seed_idx], np.arange(n_sources))


# Compute inverse solution and for each epoch. By using "return_generator=True"
# stcs will be a generator object instead of a list. This allows us so to
# compute the Coherence without having to keep all source estimates in memory.

# Now we are ready to compute the Coherence in the alpha and beta band.
# fmin and fmax specify the lower and upper freq. for each band, resp.
        fmin = (4., 8., 13., 31.)
        fmax = (7., 12., 30., 45.)
        tminw = .150
        Tn=1000*tmin
        tmaxw= .350
        Tx=1000*tmax
        sfreq = raw.info['sfreq']  # the sampling frequency

# Now we compute connectivity. To speed things up, we use 2 parallel jobs
# and use mode='multitaper', mt_adaptive=True, which uses a FFT with a Hanning window
# to compute the spectra (instead of multitaper estimation, which has a
# lower variance but is slower). By using faverage=True, we directly
# average the Coherence in the alpha and beta band, i.e., we will only
# get 2 frequency bins
	
        print "spectral connectivity coherence"
        coh, freqs1, times1, n_epochs1, n_tapers1 = spectral_connectivity(stcs, method='coh', mode='multitaper', mt_adaptive=True, indices=indices,sfreq=sfreq, fmin=fmin, fmax=fmax, faverage=True, tmin=tminw, tmax=tmaxw, n_jobs=4)
##spectral_connectivity(data, method='coh', indices=None, sfreq=6.283185307179586, mode='multitaper', fmin=None, fmax=inf, fskip=0, faverage=False, tmin=None, tmax=None, mt_bandwidth=None, mt_adaptive=False, mt_low_bias=True, cwt_frequencies=None, cwt_n_cycles=7, block_size=1000, n_jobs=1, verbose=None)
        print('coh fin ')

# Generate a SourceEstimate with the Coherence. This is simple since we
# used a single seed. For more than one seeds we would have to split coh.
# Note: We use a hack to save the frequency axis as time
        print "use a hack to save the frequency axis as time"
        tmin1 = 1#np.mean(freqs[0])
        tstep1 = 1#np.mean(freqs[1]) - tmin
        coh_stc = mne.SourceEstimate(coh, vertices=stc_coh.vertices, tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
        data_out = data_path + meg +  'firstMorphed_SemLoc_icaclean_equalized_Concrete_Coherence_ThetaAlphaBetaGamma_150_350_' + label_name[0:-3]
        print data_out
        coh_stc.save(data_out)
	
#############
#
        print "spectral connectivity plv"
        plv, freqs, times, n_epochs, n_tapers = spectral_connectivity(stcs, method='plv', mode='multitaper', mt_adaptive=True, indices=indices,sfreq=sfreq, fmin=fmin, fmax=fmax, faverage=True, tmin=tminw, tmax=tmaxw, n_jobs=4)
        #spectral_connectivity(data, method='coh', indices=None, sfreq=6.283185307179586, mode='multitaper', fmin=None, fmax=inf, fskip=0, faverage=False, tmin=None, tmax=None, mt_bandwidth=None, mt_adaptive=False, mt_low_bias=True, cwt_frequencies=None, cwt_n_cycles=7, block_size=1000, n_jobs=1, verbose=None)
        
        print('plv fin ')
        
        # Generate a SourceEstimate with the Coherence. This is simple since we
        # used a single seed. For more than one seeds we would have to split coh.
        # Note: We use a hack to save the frequency axis as time
        print "use a hack to save the frequency axis as time"
        tmin1 = 1#np.mean(freqs[0])
        tstep1 = 1#np.mean(freqs[1]) - tmin
        plv_stc = mne.SourceEstimate(plv, vertices=stc_coh.vertices, tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
        data_out = data_path + meg +  'firstMorphed_SemLoc_icaclean_equalized_Concrete_PLV_ThetaAlphaBetaGamma_150_350_' + label_name[0:-3]
        print data_out
        plv_stc.save(data_out)

###############
	print "spectral connectivity ppc"
        plv, freqs, times, n_epochs, n_tapers = spectral_connectivity(stcs, method='ppc', mode='multitaper', mt_adaptive=True, indices=indices,sfreq=sfreq, fmin=fmin, fmax=fmax, faverage=True, tmin=tminw, tmax=tmaxw, n_jobs=4)
        #spectral_connectivity(data, method='coh', indices=None, sfreq=6.283185307179586, mode='multitaper', fmin=None, fmax=inf, fskip=0, faverage=False, tmin=None, tmax=None, mt_bandwidth=None, mt_adaptive=False, mt_low_bias=True, cwt_frequencies=None, cwt_n_cycles=7, block_size=1000, n_jobs=1, verbose=None)
        
        print('ppc fin ')
        
        # Generate a SourceEstimate with the Coherence. This is simple since we
        # used a single seed. For more than one seeds we would have to split coh.
        # Note: We use a hack to save the frequency axis as time
        print "use a hack to save the frequency axis as time"
        tmin1 = 1#np.mean(freqs[0])
        tstep1 = 1#np.mean(freqs[1]) - tmin
        plv_stc = mne.SourceEstimate(plv, vertices=stc_coh.vertices, tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
        data_out = data_path + meg +  'firstMorphed_SemLoc_icaclean_equalized_Concrete_PPC_ThetaAlphaBetaGamma_150_350_' + label_name[0:-3]
        print data_out
        plv_stc.save(data_out)
