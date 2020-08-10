# Authors: Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
#          Denis A. Engemann <denis.engemann@gmail.com>
#
# License: BSD (3-clause)

print(__doc__)

import mne
import warnings
warnings.filterwarnings("ignore", category=UserWarning)
from mne import io
from mne.io import Raw
import numpy as np
import pylab as pl

###############################################################################
# Set parameters
data_path = '/imaging/rf02/TypLexMEG/'
event_path= '/group/erp/data/olaf.hauk/MEG/TypLexM/data/MF22_new/'

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
for ii, meg in enumerate(list_all):
	raw_fname = data_path + meg + 'imlex_pw_ssstf_fft_1_48_clean_ica_raw.fif'
	event_fname = event_path + meg +  'imlex_pw_raw_ssstnew.txt'
	tmin, tmax = -0.5, 0.70
	events = mne.read_events(event_fname)

	if ii < 2:
	    	word_ids=np.array([1]); #print events[events[:,2]==1,2].shape#+events[events[:,2]==2,2].shape
	    	nonword_ids=np.array([2]); #print events[events[:,2]==2,2].shape#+events[events[:,2]==6,2].shape
	    	events=mne.merge_events(events,word_ids,601,replace_events=True)
	    	events=mne.merge_events(events,nonword_ids,602,replace_events=True)
	    	event_ids = {'HIm': 601, 'LIm': 602}

	else:
	    	word_ids=np.array([7]); #print events[events[:,2]==7,2].shape#+events[events[:,2]==8,2].shape
	    	nonword_ids=np.array([8]); #print events[events[:,2]==8,2].shape#+events[events[:,2]==6,2].shape
	    	events=mne.merge_events(events,word_ids,601,replace_events=True)
	    	events=mne.merge_events(events,nonword_ids,602,replace_events=True)
	    	event_ids = {'HIm': 601, 'LIm': 602}

	#   Setup for reading the raw data
	raw = io.Raw(raw_fname)
	stim_delay=0.034 # delay in s ((NOTE! not in the example, taken from Olaf's script. Rezvan))
	events[:,0] = events[:,0] + np.round( raw.info['sfreq']*stim_delay ) #((NOTE! not in the example, taken from Olaf's script. Rezvan))

	#   Plot raw data
	#fig = raw.plot(events=events)

	#   Set up pick list: EEG + STI 014 - bad channels (modify to your needs)
	include = []  # or stim channels ['STI 014']
	exclude = []#raw.info['bads'] # bads

	# pick EEG and MEG channels
	picks = mne.pick_types(raw.info, meg=True, eeg=True, stim=False, eog=True,
		               include=include, exclude=exclude)
	# Read epochs
	if ii in np.array([3,4,13]):
         reject = dict(eeg=150e-6, grad=200e-12, mag=4e-12)
	else:
         reject = dict(eeg=120e-6, grad=200e-12, mag=4e-12)
	epochs = mne.Epochs(raw, events, event_ids, tmin, tmax, picks=picks,
		            baseline=(None, 0), reject=reject, preload=True)
	epochs_words = epochs['HIm']
	epochs_pwords = epochs['LIm']
	print ii; print epochs_words.get_data().shape[0]/float(events[events[:,2]==601,2].shape[0])
	print epochs_pwords.get_data().shape[0]/float(events[events[:,2]==602,2].shape[0])



	#epochs.plot()
	## averaging
#	out_evoked_fname = data_path + meg + 'imlex_pw_ssstf_fft_1_48_clean_ica_raw-ave.fif'
#	print "averaging"
#	evokeds = [epochs[name].average() for name in ['HIm', 'LIm']]  # average epochs and get an Evoked dataset.
#
#	#evoked.save(out_evoked_fname)  # save evoked data to disk
#	print out_evoked_fname
	#mne.write_evokeds(out_evoked_fname, evokeds)
	"""
	###############################################################################
	### View evoked response
	import matplotlib.pyplot as plt
	#plt.clf()
	#ax = plt.figure()
	##evokeds[0].plot()
	#plt.title('EEG evoked potential, abstract word trials')
	#plt.ylabel('Potential (uV)')
	#ax = plt.subplot(2, 1, 2)
	#plt.figure()
	##evokeds[1].plot()
	#plt.title('EEG evoked potential, concrete word trials')
	#plt.ylabel('Potential (uV)')
	#plt.show()

	##############
	print "TF analysis"
	# frequency bands with variable number of cycles for wavelets
	frequencies = np.arange(8, 45, 1)
	n_cycles = frequencies / float(7)
	n_cycles[frequencies<=20] = 2
	n_cycles[frequencies<=20] += np.arange(0,1,0.077)
	Fs = raw.info['sfreq']  # sampling in Hz
	from mne.time_frequency import induced_power
	epochs_data = epochs['abs_wrd'].get_data() 
	power, phase_lock = induced_power(epochs_data, Fs=Fs, frequencies=frequencies, n_cycles=n_cycles, n_jobs=1)

	print "plot power TF "
	pl.figure()
	pl.subplot(3,1,1)
	pl.imshow(power[20, :, :], aspect='auto', origin='lower')
	pl.title('Power')
	pl.colorbar
	## plot phase lock TF
	pl.subplot(3,1,2)
	pl.imshow(phase_lock[20, :, :], aspect='auto', origin='lower')
	pl.title('Phase lock')
	pl.colorbar
	pl.show()
	"""
