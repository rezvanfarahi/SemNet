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
sys.path.insert(1,'/imaging/local/software/mne_python/latest_v0.9full')
#sys.path.insert(1,'/imaging/rf02/TypLexMEG/mne_python_v8')


# for qsub
# (add ! here if needed) /imaging/local/software/anaconda/latest/x86_64/bin/python
sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
###


import numpy as np
import mne
from mne import io
from mne.minimum_norm import apply_inverse, read_inverse_operator
from mne.epochs import equalize_epoch_counts

###############################################################################
# Set parameters
main_path = '/imaging/rf02/Semnet/'
data_path = '/imaging/rf02/Semnet'	# where subdirs for MEG data are
event_path = '/imaging/rf02/Semnet'
print sys.argv
subject_inds = []
for ss in sys.argv[1:]:
 subject_inds.append( int( ss ) )

subject_inds=[0, 1,2,3,4,5,6,7,8,9,10,11,12,13]#]
print "subject_inds:"
print subject_inds
print "No rejection"

list_all =  ['/meg16_0045/160303/', #8
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


ll = []
for ss in subject_inds:
 ll.append(list_all[ss])

print "ll:"
print ll
n_subjects=14
stim_delay = 0.034 # delay in s ((NOTE! not in the example, taken from Olaf's script. Rezvan))
tmin, tmax = -0.1, 0.5
nepochs=np.zeros(25)
all_bads=np.zeros(25)
import matplotlib.pyplot as plt
for ii, meg in enumerate(ll):
    print subject_inds[ii]
    #cov_path=data_path+meg+'noise_cov'
    #noise_cov=mne.read_cov(cov_path)
    raw_fname = data_path + meg + 'block_localisers_tsss_filt_pnt1_48_ica_raw.fif'#'block_milk_tsss_filt_pnt1_48_ica_raw.fif' #clean_ssp_
    event_fname = event_path + meg + 'block_localisers_tsss_ds_raw-eve.fif'
    subjects_dir = data_path 
    
    tmin = -0.1
    tmax = 0.5  # Use a lower tmax to reduce multiple comparisons
    #   Setup for reading the raw data
    raw = io.Raw(raw_fname)#, preload=True)
    events = mne.read_events(event_fname)
    stim_delay = 0.034 # delay in s
    events[:,0]= np.round( raw.info['sfreq']*stim_delay )+events[:,0]
    ###############################################################################
    # Read epochs for all channels, removing a bad one
#    raw.plot(scalings=dict(eeg=100e-6))
    print "epochs"
    picks = mne.pick_types(raw.info, eeg=True, meg=True, eog=True, exclude='bads')
    event_id = {'audio': 1,  'colour': 3, 'grey': 4, 'shapes': 7, 'scrambled': 8}#, 'button': 4096}
    color = {1: 'blue', 3: 'red', 4: 'green', 7: 'c', 8: 'black'}#, 4096: 'yellow'}
    if subject_inds[ii] in np.array([0]):#np.array([4,8,9]):
        reject = dict(eeg=150e-6, grad=200e-12, mag=4e-12)#, eog=150e-6)
    else:
        reject = dict(eeg=120e-6, grad=200e-12, mag=4e-12)#, eog=150e-6)
    
    epochs = mne.Epochs(raw, events, event_id, tmin, tmax, picks=picks,
    		             baseline=(None,0), proj=True, reject=reject, preload=True)
    for eegcnt in range(71):
        if eegcnt<10:
            thiseegbad=sum(epochs.drop_log,[]).count('EEG00'+str(eegcnt))
            if thiseegbad>=35:
                raw.info['bads'].append(u'EEG00'+str(eegcnt))                
        else:
            thiseegbad=sum(epochs.drop_log,[]).count('EEG0'+str(eegcnt))
            if thiseegbad>=35:
                raw.info['bads'].append(u'EEG0'+str(eegcnt))
    print raw.info['bads']
    all_bads[ii]=len(raw.info['bads'])
#    if raw.info['bads']:
#         raw.interpolate_bads(reset_bads=True)
    picks = mne.pick_types(raw.info, eeg=True, meg=True, eog=True, exclude='bads')
    event_id=7# {'colour': 3, 'grey': 4, 'shapes': 7, 'scrambled': 8}
    epochs1 = mne.Epochs(raw, events, event_id, tmin, tmax, picks=picks,
    		             baseline=(None,0), proj=True, reject=reject, preload=True) 
    event_id=8# {'colour': 3, 'grey': 4, 'shapes': 7, 'scrambled': 8}
    epochs2 = mne.Epochs(raw, events, event_id, tmin, tmax, picks=picks,
    		             baseline=(None,0), proj=True, reject=reject, preload=True)                 
    print epochs
    ei=1
    event_names=['audio','colour_shape','button']
    directory='/imaging/rf02/Semnet/results/localisers/evoked/'+event_names[ei]+'/shape_sub'
    import os
    if not os.path.exists(directory):
        os.makedirs(directory)
    evoked1=epochs1.average()
    evoked2=epochs2.average() 
    evoked_sub=np.subtract(evoked1,evoked2)
    fig=evoked_sub.plot(exclude='bads', show=False)
    figname_out=directory+'/subj'+str(subject_inds[ii])+'.jpg'
    fig.savefig(figname_out)
    print epochs.events.shape[0]/100.
#    print sorted(epochs.drop_log)
#    nepochs[ii]=epochs.events.shape[0]/350.
#    path_epochs='./epochs_milk.mat'
    plt.close("all")

#import scipy.io as scio
#scio.savemat(path_epochs,{'nepochs':nepochs})
#path_bads='./bads_odour.mat'
#scio.savemat(path_bads,{'all_bads':all_bads})
#                   
#                   
#                   
#	#    Equalize trial counts to eliminate bias (which would otherwise be
#	#    introduced by the abs() performed below)
#
#	###############################################################################
#    print "Transform to source space"
#    
#    fname_inv = inv_path + meg + inv_fname
#    
#    method = "MNE"  # use dSPM method (could also be MNE or sLORETA)
#    inverse_operator = read_inverse_operator(fname_inv)
#    sample_vertices = [s['vertno'] for s in inverse_operator['src']]
#    #    Let's average and compute inverse, resampling to speed things up
#    evoked=epochs.average()    
##    epdata=evoked.data 
##    em=np.abs(epdata)
##    #emax=np.mean(em[:,:,600:900],axis=2)#
##    sa=em[:,550:950]#.max(axis=2)
##    #bm=np.sqrt(np.mean(epdata[:,:,100:500]**2,axis=2))
##    bv=np.var(em[:,100:500],axis=1)
##    edv=sa.copy()
##    for ii in range(sa.shape[1]):
##        edv[:,ii]=(sa[:,ii]**2)/bv
##
##    ediv=np.sqrt(edv)
#    snr = 3.0#np.mean(ediv)
#    print snr
#    lambda2 = 1.0 / snr ** 2
#    
#    evoked1 = epochs1.average()
#    
#    #evoked1 = whiten_evoked(evoked11, noise_cov)
#    #evoked1.resample(50)
#    condition1 = apply_inverse(evoked1, inverse_operator, lambda2, method)
#    
#    evoked2 = epochs2.average()
#    #evoked2 = whiten_evoked(evoked21, noise_cov)
#    #evoked2.resample(50)
#    condition2 = apply_inverse(evoked2, inverse_operator, lambda2, method)
#    
#    subject_from = subjects[subject_inds[0]]
#    subject_to = 'fsaverage'
#    vertices_to = [np.arange(10242), np.arange(10242)]
#    morphmat1=mne.compute_morph_matrix(subject_from, subject_to, condition1.vertices, vertices_to, subjects_dir=data_path)
#    morphmat2=mne.compute_morph_matrix(subject_from, subject_to, condition2.vertices, vertices_to, subjects_dir=data_path)
#    stc1=condition1.morph_precomputed(subject_to,vertices_to,morphmat1, subject_from)
#    stc2=condition2.morph_precomputed(subject_to,vertices_to,morphmat2, subject_from)
## 
#    data_out = data_path + meg + 'firstMorphed_ico_SemLoc_ica_Concrete_Source_Evoked_m500_700'
#    stc1.save(data_out)
#    data_out = data_path + meg + 'firstMorphed_ico_SemLoc_ica_Abstract_Source_Evoked_m500_700'
#    stc2.save(data_out)
