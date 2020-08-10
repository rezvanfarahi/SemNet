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
subject_inds = [5]
for ss in sys.argv[1:]:
   subject_inds.append( int( ss ) )

#subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14]
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
	in_fif_fname = main_path + meg + 'semloc_raw_ssst.fif'#fft48_
	out_name=data_path+meg+ 'semloc_ssstf_fft48_clean_raw.fif'
	raw = mne.io.Raw(in_fif_fname, preload=True)
	include = []
	exclude = raw.info['bads'] # bads
	picks = mne.pick_types(raw.info, meg=True, eeg=True, eog=True, stim=False, include=include, exclude=exclude)
	raw.filter(l_freq=1, h_freq=48, l_trans_bandwidth=0.1, picks=picks, method='fft')

	#ecg_events, ch_ecg, avg_pules= mne.preprocessing.find_ecg_events(raw, event_id=5999, ch_name=None, tstart=0.0, l_freq=5, h_freq=35, qrs_threshold='auto',filter_length='10s', verbose=None)

	#eog_events= mne.preprocessing.find_eog_events(raw, event_id=5998, l_freq=1, h_freq=10, filter_length='10s', ch_name=None, tstart=0,verbose=None)

	proj_eog, eog_events2=mne.preprocessing.compute_proj_eog(raw, tmin=-0.2, tmax=0.2, n_grad=2, n_mag=2, n_eeg=2, filter_length='10s', n_jobs=4, reject={'eog':10000, 'eeg': 0.0005, 'grad': 2e-10, 'mag': 3e-12}, flat=None,bads=exclude, avg_ref=False, no_proj=False, event_id=5998, eog_l_freq=1, eog_h_freq=10, tstart=0.0, filter_method='fft', iir_params=None, copy=True, verbose=None)
	
	proj_ecg, ecg_events2= mne.preprocessing.compute_proj_ecg(raw,  tmin=-0.2, tmax=0.4, n_grad=1, n_mag=2, n_eeg=2, average=False, filter_length='10s', n_jobs=4, ch_name='MEG1531', reject={'eog': 0.00025, 'eeg': 5e-05, 'grad': 2e-10, 'mag': 3e-12}, flat=None, bads=exclude, avg_ref=False, no_proj=False, event_id=5999, ecg_l_freq=5, ecg_h_freq=35, tstart=0.0,qrs_threshold='auto', filter_method='fft', iir_params=None, copy=True, verbose=None)

	raw.add_proj(proj_ecg[1:])
	raw.add_proj(proj_eog[1:])
	raw.apply_proj()
	tmin, tmax = -0.5, 0.7
	reject = dict(eeg=120e-6, grad=200e-12, mag=4e-12)
	event_ids = {'cncrt_wrd': 1, 'abs_wrd': 2}	
	#event_fname = event_path + meg + 'semloc_ssstf_fft48_raw-eve.fif'
	#print "Reading events from" + event_fname
	#events = mne.read_events(event_fname)
	events = mne.find_events(raw, stim_channel='STI101')
	stim_delay=0.034
	events[:,0] += np.round( raw.info['sfreq']*stim_delay )

## get configuration info
	
## epoching
	print "mne.Epochs()"
 
    # get epochs from raw
	epochs = mne.Epochs(raw, events, event_ids, tmin, tmax, picks=picks, baseline=(None, 0), reject=reject, preload=True)
	evokeds = [epochs[cond].average() for cond in ['cncrt_wrd', 'abs_wrd']]
	import matplotlib.pyplot as plt
	plt.clf()
	evokeds[0].plot()

	#raw.apply_proj()
	some_picks = mne.pick_types(raw.info, meg=False, eeg=True, eog=False, stim=False, include=include, exclude=exclude)
	data=raw.copy()
	data, times = data[some_picks,2000:200000]
	data=data.transpose();

	# save 150s of MEG data in FIF file
	#raw.save('sample_audvis_meg_raw.fif', tmin=0, tmax=150, picks=picks,
	 #       overwrite=True)

	###############################################################################
	# Show MEG data
	plt.plot(data)
	plt.show()
"""
	print meg
	in_path=data_path+meg
	in_fif_fname = data_path + meg + 'semloc_ssstf_fft48_raw.fif'
	raw_in = mne.io.Raw(in_fif_fname)
	prefix = in_fif_fname[:-8]	 
	out_fif_fname = prefix + '_clean_ecg_eog_raw.fif'
	ecg_proj_fname = prefix + '_ecg-proj.fif'
	eog_proj_fname = prefix + '_eog-proj.fif'
	ecg_event_fname = prefix + '_ecg-eve.fif'
	eog_event_fname = prefix + '_eog-eve.fif'
	ecg_events, _, _ = mne.preprocessing.find_ecg_events(raw_in)
	mne.write_events(ecg_event_fname, ecg_events)
	print("Writing ECG events in %s" % ecg_event_fname)
	eog_events = mne.preprocessing.find_eog_events(raw_in)
	print("Writing EOG events in %s" % eog_event_fname)
	mne.write_events(eog_event_fname, eog_events)


'mne_process_raw', '--cd', in_path, '--raw', in_fif_fname,
                   '--events', ecg_event_fname, '--makeproj',
                   '--projtmin', '-0.08', '--projtmax', '0.08',
                   '--saveprojtag', '_ecg-proj', '--projnmag', '2',
                   '--projngrad', '1', '--projevent', '999', '--highpass', '5',
                   '--lowpass', '35', '--projmagrej', '4000',
                   '--projgradrej'

print('Computing EOG projector')
command = ('mne_process_raw', '--cd', in_path, '--raw', in_fif_fname,
           '--events', eog_event_fname, '--makeproj',
           '--projtmin', '-0.15', '--projtmax', '0.15',
           '--saveprojtag', '_eog-proj', '--projnmag', '2',
           '--projngrad', '2', '--projevent', '998', '--lowpass', '35',
           '--projmagrej', '4000', '--projgradrej',


'mne_process_raw', '--cd', in_path, '--raw', in_fif_fname,
                   '--proj', in_fif_fname, '--projoff', '--save',
                   out_fif_fname, '--filteroff',
                   '--proj', ecg_proj_fname, '--proj', eog_proj_fname

"""
   
