print(__doc__)

# Authors: Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
#
# License: BSD (3-clause)

import numpy as np
import matplotlib.pyplot as plt

import mne
from mne import io


###############################################################################
# Set parameters
data_path = '/imaging/rf02/TypLexMEG/'
event_path = '/imaging/rf02/TypLexMEG/'
meg='meg11_0118/110509/'
#stim_delay = 0.034 # delay in s
raw_fname = data_path + meg + '/semloc_ssstf_fft49_clean_ecg_eog_raw.fif'
event_fname = event_path + meg + 'semloc_ssstf_fft49_raw-eve.fif'
event_id, tmin, tmax = 1, -0.5, 0.7

# Setup for reading the raw data
raw = io.Raw(raw_fname)
events = mne.read_events(event_fname)

# Set up pick list: EEG + MEG - bad channels (modify to your needs)
picks = mne.pick_types(raw.info, meg=True, eeg=True, stim=True, eog=True,
                       exclude='bads')

# Read epochs
epochs = mne.Epochs(raw, events, event_id, tmin, tmax, proj=True,
                    picks=picks, baseline=(None, 0), preload=True,
                    reject=dict(eeg=120e-6, eog=150e-6, grad=200e-12, mag=4e-12))

###############################################################################
# Show event related fields images

# and order with spectral reordering
# If you don't have scikit-learn installed set order_func to None
from sklearn.cluster.spectral import spectral_embedding
from sklearn.metrics.pairwise import rbf_kernel


def order_func(times, data):
    this_data = data[:, (times > 0.0) & (times < 0.500)]
    this_data /= np.sqrt(np.sum(this_data ** 2, axis=1))[:, np.newaxis]
    return np.argsort(spectral_embedding(rbf_kernel(this_data, gamma=1.),
                      n_components=1, random_state=0).ravel())

good_pick = 20  # channel with a clear evoked response
bad_pick = 21  # channel with no evoked response

plt.close('all')
mne.viz.plot_image_epochs(epochs, [good_pick, bad_pick], sigma=0.5, vmin=-150,
                          vmax=350, colorbar=True, order=order_func, show=True)
