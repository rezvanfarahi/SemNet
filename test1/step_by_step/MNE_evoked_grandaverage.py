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
sys.path.insert(1,'/imaging/local/software/mne_python/latest_v0.9')
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
data_path = '/imaging/rf02/TypLexMEG/'	# where subdirs for MEG data are
event_path = '/group/erp/data/olaf.hauk/MEG/TypLexM/data/MF22_new/'
#the following four lines are for thw parallel processing, will not affect the script
print sys.argv
subject_inds = [0]
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
tmin, tmax = -0.5, 0.7 #epoch length: from 500ms prestimulus to 700 peristimulus
event_ids = {'cncrt_wrd': 1, 'abs_wrd': 2}
evoked=range(n_subjects)
for ii, meg in enumerate(ll):
    print ii
    #cov_path=data_path+meg+'noise_cov'
    #noise_cov=mne.read_cov(cov_path)
    raw_fname = data_path + meg + 'semloc_ssstf_fft_1_48_clean_ica_raw.fif' #clean_ssp_
    event_fname = event_path + meg + 'semloc_raw_ssstnew.txt'
    subjects_dir = data_path 
    
    tmin = -0.5
    tmax = 0.7  
    #   Setup for reading the raw data
    raw = io.Raw(raw_fname)
    events = mne.read_events(event_fname)
    stim_delay = 0.034 # delay in s
    events[:,0] = events[:,0]+ np.round( raw.info['sfreq']*stim_delay )   
    ###############################################################################
    # Read epochs for all channels, removing a bad one
    print "epochs"
    picks = mne.pick_types(raw.info, eeg=True, meg=True, eog=False, exclude='bads')
    event_id = 1  # concrete
    if subject_inds[0]==3:
        reject = dict(eeg=150e-6, grad=200e-12, mag=4e-12) #eog=150e-6,
    else:
        reject = dict(eeg=120e-6, grad=200e-12, mag=4e-12) #eog=150e-6,
    epochs1 = mne.Epochs(raw, events, event_id, tmin, tmax, picks=picks,
    		             baseline=(None, 0), proj=True, reject=reject, preload=True)
    
    event_id = 2  # abstract
    epochs2 = mne.Epochs(raw, events, event_id, tmin, tmax, picks=picks, baseline=(None, 0), proj=True, reject=reject, preload=True)
    event_id={'cncrt_wrd': 1, 'abs_wrd': 2}
    epochs = mne.Epochs(raw, events, event_id, tmin, tmax, picks=picks,baseline=(None, 0), proj=True, reject=reject, preload=True)
    equalize_epoch_counts([epochs1, epochs2],method='mintime')#    Equalize trial counts to eliminate bias 
	
	###############################################################################
    
    evoked[ii]=epochs.average()    
evoked_ga=mne.grand_average(evoked)
out_file=data_path+'Evoked_GrandAverage'
evoked_ga.save(out_file)
   