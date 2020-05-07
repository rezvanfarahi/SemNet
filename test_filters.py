# -*- coding: utf-8 -*-
"""
Created on Sat Sep 16 14:36:13 2017

@author: rf02
"""
print __doc__

import sys
#sys.path.insert(1,'/imaging/local/software/mne_python/v0.11')
sys.path.insert(1,'/imaging/local/software/mne_python/latest_v0.15')

import numpy as np
import mne

import numpy as np
from scipy import signal, fftpack
import matplotlib.pyplot as plt

from mne.time_frequency.tfr import morlet
from mne.viz import plot_filter, plot_ideal_filter


sfreq = 1000.
f_p = 40.
flim = (1., sfreq / 2.)  # limits for plotting
nyq = sfreq / 2.  # the Nyquist frequency is half our sample rate
freq = [0, f_p, f_p, nyq]
gain = [1, 1, 0, 0]

dur = 10.
center = 2.
morlet_freq = f_p
tlim = [center - 0.2, center + 0.2]
tticks = [tlim[0], center, tlim[1]]
flim = [20, 70]

x = np.zeros(int(sfreq * dur) + 1)
blip = morlet(sfreq, [morlet_freq], n_cycles=7)[0].imag / 20.
n_onset = int(center * sfreq) - len(blip) // 2
x[n_onset:n_onset + len(blip)] += blip
x_orig = x.copy()

rng = np.random.RandomState(0)
x += rng.randn(len(x)) / 1000.
x += np.sin(2. * np.pi * 60. * np.arange(len(x)) / sfreq) / 2000.

transition_band = 0.25 * f_p
f_s = f_p + transition_band
filter_dur = 6.6 / transition_band  # sec
n = int(sfreq * filter_dur)
freq = [0., f_p, f_s, sfreq / 2.]
gain = [1., 1., 0., 0.]
# This would be equivalent:
# h = signal.firwin2(n, freq, gain, nyq=sfreq / 2.)
h = mne.filter.create_filter(x, sfreq, l_freq=None, h_freq=f_p)
x_shallow = np.convolve(h, x)[len(h) // 2:]

plot_filter(h, sfreq, freq, gain, 'MNE-Python 0.14 default', flim=flim)