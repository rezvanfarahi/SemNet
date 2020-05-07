# -*- coding: utf-8 -*-
"""
Created on Wed Apr  6 11:06:01 2016

@author: rf02
"""

"""
Preprocessing MEG and EEG data using filters and ICA
====================================================

This examples filters MEG and EEG data and subsequently computes separate
ICA solutions for MEG and EEG. An html report is created, which includes
diagnostic plots of the preprocessing steps.
"""
# Author: Denis A. Engemann <denis.engemann@gmail.com>
# License: BSD (3-clause)

import os.path as op
import sys
sys.path.insert(1,'/imaging/local/software/mne_python/latest_v0.11')

import numpy as np
import mne


from meeg_preprocessing import check_apply_filter, compute_ica
from meeg_preprocessing.utils import get_data_picks, setup_provenance

from mne import io
from config import (
    fname,
    subject,
    n_jobs,
    filter_params,
    notch_filter_params,
    plot_fmin,
    plot_fmax,
    n_components,
    ica_reject,
    n_max_ecg,
    n_max_eog,
    decim,
    ica_meg_combined,
    img_scale,
)

# configure logging + provenance tracking magic
report, run_id, results_dir, logger = setup_provenance(
    script=__file__, results_dir='results')

raw = io.Raw(fname, preload=True)

include = []  # or stim channels ['STI 014']
raw.info['bads'] += ['EEG 053']  # bads + 1 more
epochs2 = mne.Epochs(raw, events, event_id, tmin, tmax, picks=picks, baseline=(None, 0), proj=True, reject=reject, preload=True)
fig, report = check_apply_filter(raw, subject=subject,
                                 filter_params=filter_params,
                                 notch_filter_params=notch_filter_params,
                                 plot_fmin=plot_fmin, plot_fmax=plot_fmax,
                                 n_jobs=n_jobs, img_scale=img_scale)

# get picks and iterate over channels
artifact_stats = dict()
for picks, ch_type in get_data_picks(raw, meg_combined=ica_meg_combined):
    ica, _ = compute_ica(raw, picks=picks,
                         subject=subject, n_components=n_components,
                         n_max_ecg=n_max_ecg, n_max_eog=n_max_eog,
                         reject=ica_reject,
                         random_state=42,
                         artifact_stats=artifact_stats,
                         decim=decim, report=report, img_scale=img_scale)
    ica.save('{}-ica.fif'.format(ch_type))

for k, v in artifact_stats.items():
    print(k, v)

report.save(  # save in automatically generated folder
    op.join(results_dir, run_id,
            'preprocessing-report-{}.html'.format(subject)),
    open_browser=True, overwrite=True)
