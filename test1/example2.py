import mne #
from mne.datasets import sample #
data_path = sample.data_path() #
raw_fname = data_path + '/MEG/sample/sample_audvis_filt-0-40_raw.fif' #
raw = mne.io.Raw(raw_fname) # 
start, stop = raw.time_as_index([100, 115])  # 100 s to 115 s data segment  #
picks = mne.pick_types(raw.info, meg=True, eeg=False, stim=True, exclude='bads') #
events = mne.find_events(raw, stim_channel='STI 014') #
event_id = dict(aud_l=1, aud_r=2)  # event trigger and conditions #
tmin = -0.2  # start of each epoch (200ms before the trigger) #
tmax = 0.5  # end of each epoch (500ms after the trigger) #
picks = mne.pick_types(raw.info, meg=True, eeg=True, eog=True, stim=False, exclude='bads') #
baseline = (None, 0)  # means from the first instant to t = 0 #
reject = dict(grad=4000e-13, mag=4e-12, eog=150e-6) #
epochs2 = mne.Epochs(raw, events, event_id, tmin, tmax, proj=True, picks=picks, baseline=baseline, preload=False, reject=reject) #
evoked = epochs2['aud_l'].average() 
#print(evoked)
#evoked.plot()

print "TF analysis"
import numpy as np
n_cycles = 2  # number of cycles in Morlet wavelet
frequencies = np.arange(7, 30, 3)  # frequencies of interest
Fs = raw.info['sfreq']  # sampling in Hz
from mne.time_frequency import induced_power
epochs_data = epochs2['aud_l'].get_data() 
power, phase_lock = induced_power(epochs_data, Fs=Fs, frequencies=frequencies, n_cycles=2, n_jobs=1)

print "inverse modeling"
from mne.minimum_norm import apply_inverse, read_inverse_operator
fname_inv = data_path + '/MEG/sample/sample_audvis-meg-oct-6-meg-inv.fif'
inverse_operator = read_inverse_operator(fname_inv) 
snr = 3.0
lambda2 = 1.0 / snr ** 2
method = "dSPM"
stc = apply_inverse(evoked, inverse_operator, lambda2, method)
fname_label = data_path + '/MEG/sample/labels/Aud-lh.label'
label = mne.read_label(fname_label)
from mne.minimum_norm import apply_inverse_raw
start, stop = raw.time_as_index([0, 15])  # read the first 15s of data
stc = apply_inverse_raw(raw, inverse_operator, lambda2, method, label, start, stop)


