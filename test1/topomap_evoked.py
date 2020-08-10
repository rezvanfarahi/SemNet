print(__doc__)

import numpy as np
import matplotlib.pyplot as plt
from mne.datasets import sample
from mne import read_evokeds

path = '/home/rf02/rezvan/TypLexMEG'
meg = '/meg10_0378/101209/'
fname = path + meg + 'semloc_abstract&concrete_evoked-ave.fif'

# load evoked and subtract baseline
condition = 'abs_wrd'
evoked = read_evokeds(fname, condition=condition, baseline=(None, 0))

# set time instants in seconds (from 50 to 150ms in a step of 10ms)
times = np.arange(0.05, 0.15, 0.01)
# If times is set to None only 10 regularly spaced topographies will be shown

# plot magnetometer data as topomaps
evoked.plot_topomap(times, ch_type='mag')

# plot gradiometer data (plots the RMS for each pair of gradiometers)
evoked.plot_topomap(times, ch_type='grad')

# plot magnetometer data as topomap at 1 time point : 100ms
# and add channel labels and title
evoked.plot_topomap(0.1, ch_type='mag', show_names=True, colorbar=False,
                    size=6, res=128, title='Abstract Words')
plt.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.88)
