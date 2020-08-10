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

print(__doc__)
# Russell's addition
import sys
sys.path.append('/imaging/local/software/python_packages/nibabel/1.3.0')
sys.path.append('/imaging/local/software/python_packages/pysurfer/v0.4')
# End

# for qsub
# (add ! here if needed) /imaging/local/software/anaconda/latest/x86_64/bin/python
sys.path.insert(1,'/imaging/local/software/mne_python/latest_v0.9full')
#sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
###


import numpy as np

import mne
from mne import io

from mne.minimum_norm import apply_inverse_epochs, read_inverse_operator


###############################################################################
# Set parameters
data_path = '/imaging/rf02/TypLexMEG/'	# where subdirs for MEG data are
inv_path = '/imaging/rf02/TypLexMEG/'
inv_fname = 'typlex_InvOp_ico5diagreg_fft_1_48_clean_ica_EMEG-inv.fif'#'InverseOperator_fft_1_48_clean_ica_EMEG-inv.fif'
event_path = '/group/erp/data/olaf.hauk/MEG/TypLexM/data/MF22_new/'
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
n_subjects=17
stim_delay = 0.034 # delay in s ((NOTE! not in the example, taken from Olaf's script. Rezvan))
tmin, tmax = -0.5, 0.7
for ii, meg in enumerate(ll):
    subject_no= subject_inds[0]; print subject_no
    raw_fname = data_path + meg + 'typlex_pw_ssstf_fft_ica_beta_raw.fif' #clean_ssp_
    event_fname = event_path + meg + 'typlex_pw_raw_ssstnew.txt'
    subjects_dir = data_path 
    
    tmin = -0.5
    tmax = 0.7  # Use a lower tmax to reduce multiple comparisons
    #   Setup for reading the raw data
    raw = io.Raw(raw_fname)
    events = mne.read_events(event_fname)
    stim_delay = 0.034 # delay in s
    events[:,0] += np.round( raw.info['sfreq']*stim_delay )
    if subject_no < 3:
	    	word_ids=np.array([2,1]); print events[events[:,2]==1,2].shape+events[events[:,2]==2,2].shape
	    	nonword_ids=np.array([6,5]); print events[events[:,2]==5,2].shape+events[events[:,2]==6,2].shape
	    	events=mne.merge_events(events,word_ids,601,replace_events=True)
	    	events=mne.merge_events(events,nonword_ids,602,replace_events=True)
	    	event_id = {'Words': 601, 'PWords': 602}

    else:
	    	word_ids=np.array([8,7]); print events[events[:,2]==7,2].shape+events[events[:,2]==8,2].shape
	    	nonword_ids=np.array([6,5]); print events[events[:,2]==5,2].shape+events[events[:,2]==6,2].shape
	    	events=mne.merge_events(events,word_ids,601,replace_events=True)
	    	events=mne.merge_events(events,nonword_ids,602,replace_events=True)
	    	event_id = {'Words': 601, 'PWords': 602}
    ###############################################################################
    # Read epochs for all channels, removing a bad one
    print "epochs"
    picks = mne.pick_types(raw.info, eeg=True, meg=True, eog=True, exclude='bads')
    event_id = 601  # word
    if subject_no in np.array([3,4,13]):
        reject = dict(eeg=150e-6, grad=200e-12, mag=4e-12) #eog=150e-6,
    else:
        reject = dict(eeg=120e-6, grad=200e-12, mag=4e-12) #eog=150e-6,
    #reject = dict(grad=200e-12, mag=4e-12)
    epochs1 = mne.Epochs(raw, events, event_id, tmin, tmax, picks=picks, baseline=(None, 0), proj=True, reject=reject, preload=True)

    event_id = 602  # nonword
    epochs2 = mne.Epochs(raw, events, event_id, tmin, tmax, picks=picks, baseline=(None, 0), proj=True, reject=reject, preload=True)
    epochs = mne.Epochs(raw, events, event_id, tmin, tmax, picks=picks,baseline=(None, 0), proj=True, reject=reject, preload=True)
               
    


	#    Equalize trial counts to eliminate bias (which would otherwise be
	#    introduced by the abs() performed below)
	#equalize_epoch_counts([epochs1, epochs2])

	###############################################################################
    print "Transform to source space"

    fname_inv = inv_path + meg + inv_fname
    epdata=epochs.get_data()
    em=np.abs(epdata)
    #emax=np.mean(em[:,:,600:900],axis=2)#
    sa=em[:,:,550:950]#.max(axis=2)
    #bm=np.sqrt(np.mean(epdata[:,:,100:500]**2,axis=2))
    bv=np.var(em[:,:,100:500],axis=2)
    edv=sa.copy()
    for ii in range(sa.shape[2]):
        edv[:,:,ii]=(sa[:,:,ii]**2)/bv
    
    ediv=np.mean(np.sqrt(edv),axis=2)
    print ediv.mean()
    snr = ediv.mean()
    lambda2 = 1.0 / snr ** 2
    method = "MNE"  # use dSPM method (could also be MNE or sLORETA)
    inverse_operator = read_inverse_operator(fname_inv)
    sample_vertices = [s['vertno'] for s in inverse_operator['src']]
    
    #    Let's average and compute inverse, resampling to speed things up
    
    evoked1 = epochs1.average()
    #evoked1.resample(50)
    condition1 = apply_inverse_epochs(epochs1, inverse_operator, lambda2, method, pick_ori="normal")
    print "cond1 fin"
    evoked2 = epochs2.average()
    #evoked2.resample(50)
    condition2 = apply_inverse_epochs(epochs2, inverse_operator, lambda2, method, pick_ori="normal")
    subject_from = subjects[subject_no]
    subject_to = 'fsaverage'
    print "cond2 fin"
    vertices_to = [np.arange(10242), np.arange(10242)]
    morphmat1=mne.compute_morph_matrix(subject_from, subject_to, condition1[0].vertices, vertices_to, subjects_dir=data_path)
    morphmat2=mne.compute_morph_matrix(subject_from, subject_to, condition2[0].vertices, vertices_to, subjects_dir=data_path)
    for ccc, con1 in enumerate(condition1):
        print ccc		
        stcf1=con1.morph_precomputed(subject_to,vertices_to,morphmat1, subject_from)
        data_out1 = data_path + meg + 'Morphed_ico_TypLex_epochs' + str(ccc) + '_icaclean_word_Source_Beta_m500_700'
        stcf1.save(data_out1)
    for ccc, con2 in enumerate(condition2):
        print ccc		
        stcf2=con2.morph_precomputed(subject_to,vertices_to,morphmat2, subject_from)
        data_out2 = data_path + meg + 'Morphed_ico_TypLex_epochs' + str(ccc) + '_icaclean_nonword_Source_Beta_m500_700'
        stcf2.save(data_out2)
	
	

