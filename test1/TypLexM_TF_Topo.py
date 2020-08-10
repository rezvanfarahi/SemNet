print(__doc__)

# Authors: Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
#          Denis Engemann <denis.engemann@gmail.com>
#
# License: BSD (3-clause)

import numpy as np
import mne
from mne import io
from mne.time_frequency import tfr_morlet
from mne.datasets import somato

###############################################################################
# Set parameters
data_path = '/home/rf02/rezvan/TypLexMEG'
meg = '/meg10_0378/101209/'
raw_fname = data_path + meg + 'semloc_raw_ssstf_raw.fif'
event_id = {'abs_wrd': 1, 'cncrt_wrd': 2}
tmin, tmax = -0.1, 0.5

# Setup for reading the raw data
raw = io.Raw(raw_fname)
baseline = (None, 0)
##events
event_path = '/home/rf02/rezvan/TypLexMEG/'    # where event files are
event_fname = event_path + meg + 'semloc_raw_ssstf-eve.fif'
print "Reading events from" + event_fname
events = mne.read_events(event_fname)
stim_delay = 0.034 # delay in s
events[:,0] += np.round( raw.info['sfreq']*stim_delay )

# picks channels
picks = mne.pick_types(raw.info, meg=True, eeg=True, eog=True, stim=False)
reject = dict(eeg=120e-6, eog=150e-6, grad=200e-12, mag=4e-12)
epochs = mne.Epochs(raw, events, event_id, tmin, tmax, picks=picks,
                    baseline=baseline, reject=reject)

###############################################################################
# Calculate power and intertrial coherence

freqs = np.arange(8, 45, 1)  # define frequencies of interest
n_cycles = freqs / float(7)
n_cycles[freqs<=20] = 2
n_cycles[freqs<=20] += np.arange(0,1,0.077)  # different number of cycle per frequency
epochs=epochs['abs_wrd']
power, itc = tfr_morlet(epochs, freqs=freqs, n_cycles=n_cycles, use_fft=False,
                        return_itc=True, decim=3, n_jobs=1)

# Baseline correction can be applied to power or done in plots
# To illustrate the baseline correction in plots the next line is commented
# power.apply_baseline(baseline=(-0.5, 0), mode='logratio')

# Inspect power
power.plot_topo(baseline=(-0.5, 0), mode='logratio', title='Average power')
power.plot([82], baseline=(-0.5, 0), mode='logratio')

import matplotlib.pyplot as plt
fig, axis = plt.subplots(1, 2, figsize=(7, 4))
power.plot_topomap(ch_type='grad', tmin=0.5, tmax=1.5, fmin=8, fmax=12,
                   baseline=(-0.5, 0), mode='logratio', axes=axis[0],
                   title='Alpha', vmin=-0.45, vmax=0.45)
power.plot_topomap(ch_type='grad', tmin=0.5, tmax=1.5, fmin=13, fmax=25,
                   baseline=(-0.5, 0), mode='logratio', axes=axis[1],
                   title='Beta', vmin=-0.45, vmax=0.45)
mne.viz.tight_layout()

# Inspect ITC
itc.plot_topo(title='Inter-Trial coherence', vmin=0., vmax=1., cmap='Reds')
