#!!!!!! check EEG channels
"""
cd /imaging/rf02/TypLexMEG/meg11_0101/110411/
mne_check_eeg_locations --file typlex_pw_ssstf_fft49_raw.fif

"""

print(__doc__)

#import mne
#from mne.datasets import sample
#from mne.io import Raw
#data_path = sample.data_path()
#fname = data_path + '/MEG/sample/sample_audvis_raw.fif'

import mne
from mne import io
from mne.io import Raw
data_path = '/home/rf02/rezvan/TypLexMEG'
meg = '/meg10_0378/101209/'
fname = data_path + meg + 'semloc_ssstf_fft49_clean_ecg_eog_raw.fif'
#out_path = fname[1:-8] + '-eve.fif'

# Reading events
raw = Raw(fname)

events = mne.find_events(raw, stim_channel='STI101')

# Writing events
out_name=data_path + meg + 'semloc_ssstf_fft49_raw-eve.fif'
#cd out_path
mne.write_events(out_name, events)

for ind, before, after in events[:5]:
    print("At sample %d stim channel went from %d to %d"
          % (ind, before, after))

# Plot the events to get an idea of the paradigm
# Specify colors and an event_id dictionary for the legend.
event_id = {'cncrt_wrd': 1, 'abs_wrd': 2,}
color = {1: 'green', 2: 'red'}
mne.viz.plot_events(events, raw.info['sfreq'], raw.first_samp, color=color, event_id=event_id)

# Plot the events to get an idea of the paradigm
# Specify colors and an event_id dictionary for the legend.
#event_id = {'abs_wrd': 1, 'cncrt_wrd': 2, 'vis_l': 3, 'vis_r': 4, 'smiley': 5,
#            'button': 32}
#color = {1: 'green', 2: 'yellow', 3: 'red', 4: 'c', 5: 'black', 32: 'blue'}

#mne.viz.plot_events(events, raw.info['sfreq'], raw.first_samp, color=color,
#                    event_id=event_id)

