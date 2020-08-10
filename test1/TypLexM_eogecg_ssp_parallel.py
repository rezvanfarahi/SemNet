print(__doc__)
print "hi! this script is for spectral functional connectivity. Very nice maps are fruits which worth waiting for I bet!"
import sys

sys.path.append('/imaging/local/software/python_packages/nibabel/1.3.0')
sys.path.append('/imaging/local/software/python_packages/pysurfer/v0.4')
sys.path.append('/imaging/local/software/python_packages/scikit-learn/v0.14.1')
#sys.path.append('/imaging/local/software/mne_python/latest')

# for qsub
# (add ! here if needed) 
sys.path.append('/imaging/local/software/EPD/latest/x86_64/bin/python')
sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
###

import os
import numpy as np
import mne
from mne import io
from mne.io import Raw
import operator
from mne.minimum_norm import (read_inverse_operator, make_inverse_operator, write_inverse_operator, apply_inverse, apply_inverse_epochs)
from mne.minimum_norm.inverse import (prepare_inverse_operator)
import sklearn
import scipy.io
from mne import find_events
from mne.connectivity import seed_target_indices, spectral_connectivity
data_path = '/imaging/rf02/TypLexMEG/' # root directory for your MEG data
os.chdir(data_path)
event_path = '/imaging/rf02/TypLexMEG/'    # where event files are
label_path = '/group/erp/data/olaf.hauk/MEG/TypLexM/labels/createdlabels_SL/'
inv_path = '/imaging/rf02/TypLexMEG/'
main_path = '/group/erp/data/olaf.hauk/MEG/TypLexM/data/MF22_new/'
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
ll = []
for ss in subject_inds:
 ll.append(list_all[ss])

print "ll:"
print ll

for meg in ll:
	print meg
	in_path=data_path+meg
	in_fif_fname = data_path + meg + 'semloc_ssstf_fft_1_48_raw.fif'#'semloc_raw_ssst.fif'
	out_name=data_path+meg+ 'semloc_ssstf_fft_1_48_clean_ssp_raw.fif'
	#out_filt=data_path+meg+ 'semloc_ssstf_fft_1_48_raw.fif'
	print "read raw"
	raw = mne.io.Raw(in_fif_fname, preload=True)
	include = []
	exclude = raw.info['bads'] # bads
	picks = mne.pick_types(raw.info, meg=True, eeg=True, eog=True, stim=False, include=include, exclude=exclude)
	#print "filtering"
	#raw.filter(l_freq=1, h_freq=48, l_trans_bandwidth=0.1, picks=picks, method='fft')
	#raw.save(out_filt, overwrite=True)

	#ecg_events, ch_ecg, avg_pules= mne.preprocessing.find_ecg_events(raw, event_id=5999, ch_name=None, tstart=0.0, l_freq=5, h_freq=35, qrs_threshold='auto',filter_length='10s', verbose=None)

	#eog_events= mne.preprocessing.find_eog_events(raw, event_id=5998, l_freq=1, h_freq=10, filter_length='10s', ch_name=None, tstart=0,verbose=None)
	print "project eog"
	proj_eog, eog_events2=mne.preprocessing.compute_proj_eog(raw, tmin=-0.2, tmax=0.2, n_grad=2, n_mag=2, n_eeg=2, filter_length='10s', n_jobs=4, reject={'eog':10000, 'eeg': 0.0005, 'grad': 2e-10, 'mag': 3e-12}, flat=None,bads=exclude, avg_ref=False, no_proj=True, event_id=5998, eog_l_freq=1, eog_h_freq=10, tstart=0.0, filter_method='fft', copy=True, verbose=None)
	print "project ecg"
	proj_ecg, ecg_events2= mne.preprocessing.compute_proj_ecg(raw,  tmin=-0.2, tmax=0.4, n_grad=2, n_mag=2, n_eeg=2, average=False, filter_length='10s', n_jobs=4, ch_name='MEG1531', reject={'eog': 0.00025, 'eeg': 5e-05, 'grad': 2e-10, 'mag': 3e-12}, flat=None, avg_ref=False, bads=exclude, no_proj=True, event_id=5999, ecg_l_freq=5, ecg_h_freq=35, tstart=0.0,qrs_threshold='auto', filter_method='fft', verbose=None)
	print "saving output"
	if proj_ecg != None:
	   raw.add_proj(proj_ecg)
	if proj_eog != None:
	   raw.add_proj(proj_eog)
	raw.save(out_name, overwrite=True)
	
