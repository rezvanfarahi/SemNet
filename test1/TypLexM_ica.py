# Authors: Denis Engemann <denis.engemann@gmail.com>
#
# License: BSD 3-clause

""" Run complete ICA for MEG and EEG

This tutorial demonstrates how to perform an entire
preprocessing workflow for one subject and for different sensor types.

1) Filtering
2) ICA for MEG
3) ICA for EEG
"""
print __doc__
# Russell added 22/02/2013
import sys
sys.path.append('/imaging/local/software/python_packages/nibabel/1.3.0')
sys.path.append('/imaging/local/software/python_packages/pysurfer/v0.3.1')
sys.path.insert(1, '/imaging/local/software/python_packages/scikit-learn/v0.14.1')
# END

# for qsub
# (add ! here if needed) /imaging/local/software/anaconda/latest/x86_64/bin/python
#sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
###


import os.path as op

import numpy as np
import matplotlib.pyplot as plt

import mne
from mne.datasets import sample
from mne.preprocessing import ICA
from mne import pick_types
from mne.preprocessing import create_eog_epochs
from mne.preprocessing import create_ecg_epochs
from mne import fiff

mne.set_log_level('INFO')

###############################################################################
# Parameters

data_path = '/imaging/rf02/TypLexMEG/meg11_0050/110307/'
raw_fname = data_path + 'semloc_ssstf_fft49_raw.fif'
raw = fiff.Raw(raw_fname, preload=True)


has_eeg = 'eeg' in raw
n_jobs = 4
dpi = 150
subject = 'avgsubject'
results_dir = '.'
filter_method = 'iir'
max_ecg = 3
max_eog = 2
ica_fname = '{ch_type}-ica.fif'
sss_filtered = True
n_components_default = 0.99
h_pass_freq = 1
l_pass_freq = 48


def save_fig(name):
    fig = plt.gcf()  # get current figure
    fig.savefig(op.join(results_dir, name), dpi=dpi)

###############################################################################
# Filtering & Checking

fig, axes = plt.subplots(1, 3 if has_eeg else 2, sharey=True, sharex=True)

# ax1, ax2, ax3, ax4 = axes.flatten()

if has_eeg:
    ax1, ax2, ax3 = axes.flatten()
else:
    ax1, ax2 = axes.flatten()

picks_mag = pick_types(raw.info, meg='mag', eeg=False)
picks_grad = pick_types(raw.info, meg='grad', eeg=False)
if has_eeg:
    picks_eeg = pick_types(raw.info, meg=False, eeg=True)


raw.plot_psds(fmin=h_pass_freq, fmax=l_pass_freq + 20, ax=ax1,
              picks=picks_mag, color='black')
line1_1 = ax1.get_lines()[0]
line1_1.set_label('mag - raw')
ax1.set_ylabel('Power (dB)')

raw.plot_psds(fmin=h_pass_freq, fmax=l_pass_freq + 20, ax=ax2,
              picks=picks_grad, color='black')
ax2.get_lines()[0].set_label('grad - raw')
ax1.set_xlabel('Frequency (Hz)')
ax1.set_ylabel('Power (dB)')

if has_eeg:
    raw.plot_psds(fmin=h_pass_freq, fmax=l_pass_freq + 20,
                  ax=ax3, picks=picks_eeg, color='black')
    ax3.get_lines()[0].set_label('eeg - raw')
    ax1.set_ylabel('Power (dB)')

raw.notch_filter(np.arange(50, 251, 50), n_jobs=n_jobs)
raw.filter(l_freq=h_pass_freq, h_freq=l_pass_freq,
           n_jobs=n_jobs, method=filter_method)

raw.plot_psds(fmin=h_pass_freq, fmax=l_pass_freq + 20,
              ax=ax1, picks=picks_mag, color='red')
ax1.get_lines()[1].set_label('mag - filt')

raw.plot_psds(fmin=h_pass_freq, fmax=l_pass_freq + 20,
              ax=ax2, picks=picks_grad, color='red')
ax2.get_lines()[1].set_label('grad - filt')

if has_eeg:
    raw.plot_psds(fmin=h_pass_freq, fmax=l_pass_freq + 20,
                  ax=ax3, picks=picks_eeg, color='red')
    ax2.get_lines()[1].set_label('eeg - filt')

fig.suptitle('Multitaper PSD')
plt.legend(loc='best')
save_fig('{}_psd_spectra.png'.format(subject))

# create ECG epochs to improve detection by correlation
picks_artifact = mne.pick_types(raw.info, meg=True, eeg=True, ecg=True,
                                eog=True)
ecg_epochs = create_ecg_epochs(raw, picks=picks_artifact)
eog_epochs = create_eog_epochs(raw, picks=picks_artifact)

for ch_type in ('eeg', 'meg'):
    if ch_type not in raw:
        continue
    if ch_type == 'meg':
        picks = mne.pick_types(raw.info, meg=True, eeg=False)
        reject = dict(mag=4e-12, grad=4000e-13)
    else:
        picks = mne.pick_types(raw.info, meg=False, eeg=True)
        reject = dict(eeg=150e-6)

    channels_picked = [raw.ch_names[k] for k in picks]
    ecg_evoked = ecg_epochs.pick_channels(channels_picked, copy=True).average()
    eog_evoked = eog_epochs.pick_channels(channels_picked, copy=True).average()
    ch_type_ = 'eeg' if ch_type == 'eeg' else 'mag'

    ######################################################################
    # Setup ICA seed decompose data, then access and plot sources.

    rank = raw.estimate_rank(picks=picks)
    if sss_filtered:
        n_components = rank
    else:  # always use rank estimate for EEG
        n_components = n_components_default if ch_type != 'eeg' else rank

    ica = mne.preprocessing.ica.ICA(n_components=n_components, max_pca_components=None)

    ##################################
    # Fit ICA model and identify bad sources

    ica.fit(raw, picks=picks, decim=3, reject=reject)

    # ECG

    ecg_inds, scores = ica.find_bads_ecg(ecg_epochs)  # inds sorted!

    if np.any(ecg_inds):
        ecg_inds = ecg_inds[:max_ecg]
        ica.plot_scores(scores, exclude=ecg_inds)  # inspect metrics used
        save_fig('{}_ica_scores_ecg_{}.png'.format(subject, ch_type))

        # indices of top five scores
        show_picks = np.abs(scores).argsort()[::-1][:5]

        # detected artifacts drawn in red (via exclude)
        ica.plot_sources(raw, show_picks, exclude=ecg_inds, start=0., stop=3.0)
        save_fig('{}_ica_sources_ecg_{}.png'.format(subject, ch_type))
        # show component sensitivites
        ica.plot_components(ecg_inds, ch_type=ch_type_, colorbar=False,
                            title='ECG')
        save_fig('{}_ica_components_ecg_{}.png'.format(subject, ch_type))

        # latent ECG sources + selection
        ica.plot_sources(ecg_evoked, exclude=ecg_inds)
        save_fig('{}_ica_overlay_ecg_evoked_{}.png'.format(subject, ch_type))

        # overlay raw and clean ECG artifacts
        ica.plot_overlay(ecg_evoked, exclude=ecg_inds)
        save_fig('{}_ica_sources_ecg_evoked_{}.png'.format(subject, ch_type))

        ica.exclude += ecg_inds  # mark first for exclusion

    # EOG

    eog_inds, scores = ica.find_bads_eog(eog_epochs)  # inds sorted!
    if np.any(eog_inds):
        eog_inds = eog_inds[:max_eog]

        ica.plot_scores(scores, exclude=eog_inds)  # inspect metrics used
        save_fig('{}_ica_scores_eog_{}.png'.format(subject, ch_type))

        # indices of top five scores
        show_picks = np.abs(scores).argsort()[::-1][:5]

        # detected artifacts drawn in red (via exclude)
        ica.plot_sources(raw, show_picks,
                         exclude=eog_inds, start=0., stop=3.0)
        save_fig('{}_ica_sources_eog_{}.png'.format(subject, ch_type))

        # show component sensitivites
        ica.plot_components(eog_inds, ch_type=ch_type_, colorbar=False,
                            title='EOG')
        save_fig('{}_ica_components_eog_{}.png'.format(subject, ch_type))

        ica.plot_sources(eog_evoked, exclude=eog_inds)
        # latent EOG sources + selection
        save_fig('{}_ica_overlay_eog_evoked_{}.png'.format(subject, ch_type))

        ica.plot_overlay(eog_evoked, exclude=eog_inds)
        # overlay raw and clean EOG artifacts
        save_fig('{}_ica_sources_eog_evoked_{}.png'.format(subject, ch_type))

        ica.exclude += eog_inds  # mark first for exclusion

    if np.any(ica.exclude):
        ica.plot_overlay(raw, start=20., stop=25.)  # EOG/ECG artifacts remain
        save_fig('{}_ica_overlay_raw_evoked_{}.png'.format(subject, ch_type))

    ica.save(op.join(results_dir, ica_fname.format(ch_type=ch_type)))

