"""

.. _tut_artifacts_correct_ica:

Artifact Correction with ICA
============================

ICA finds directions in the feature space
corresponding to projections with high non-Gaussianity. We thus obtain
a decomposition into independent components, and the artifact's contribution
is localized in only a small number of components.
These components have to be correctly identified and removed.

If EOG or ECG recordings are available, they can be used in ICA to
automatically select the corresponding artifact components from the
decomposition. To do so, you have to first build an Epoch object around
blink or heartbeat event.
"""
import sys
sys.path.insert(1,'/imaging/local/software/mne_python/v0.11')
sys.path.insert(1,'/home/rf02/semnet/meeg-preprocessing')
sys.path.append('/home/rf02/rezvan/test1')
#sys.path.append('/imaging/local/software/python_packages/scikit-learn/v0.14.1')
#sys.path.insert(1,'/imaging/local/software/python_packages/scikit-learn/v0.16.1')
sys.path.insert(1,'/imaging/rf02/scikit-learn-0.16.1')
sys.path.insert(1,'/home/rf02/.local/lib/python2.7/site-packages')
# for qsub
# (add ! here if needed) 
#sys.path.append('/imaging/local/software/anaconda/latest/x86_64/bin/python')
#sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
###
import os.path as op
import mne
reload(mne)
import sklearn
#import pylab as pl
#from mne import (fiff, read_evokeds, equalize_channels,read_proj, read_selection)

#from mne import filter
#import meeg_preprocessing
from preprocessing3 import compute_ica#from meeg_preprocessing.preprocessing import compute_ica

from utils import get_data_picks, setup_provenance
from mne.preprocessing import read_ica
import numpy as np


from mne.preprocessing import ICA
from mne.preprocessing import create_eog_epochs, create_ecg_epochs

data_path = '/imaging/rf02/Semnet/'
print sys.argv
subject_inds = []
for ss in sys.argv[1:]:
   subject_inds.append( int( ss ) )

subject_inds=[3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]#0, 
print "subject_inds:"
print subject_inds
print "No rejection"

list_all =  ['/meg16_0030/160216/', #0
            '/meg16_0032/160218/', #1
            '/meg16_0033/160218/', #2
            '/meg16_0034/160219/', #3
            '/meg16_0035/160222/', #4
            '/meg16_0039/160225/', #5
            '/meg16_0041/160226/', #6
            '/meg16_0042/160229/', #7
            '/meg16_0045/160303/', #8
            '/meg16_0047/160304/', #9
            '/meg16_0052/160310/', #10
            '/meg16_0056/160314/',#11
            '/meg16_0069/160405/',#12 
            '/meg16_0070/160407/', #13
            '/meg16_0071/160407/', #14
            '/meg16_0072/160408/', #15
            '/meg16_0073/160411/', #16
            '/meg16_0075/160411/', #17
            '/meg16_0078/160414/', #18
            '/meg16_0082/160418/', #19
            '/meg16_0086/160422/', #20
            '/meg16_0097/160512/', #21 
            '/meg16_0122/160707/', #22 LD
            '/meg16_0123/160708/', #23 LD
            '/meg16_0125/160712/', #24 LD
            ]

subjects=['S000',
'S001',
'S002',
'S003',
'S004',
'S005',
'S006',
'S007',
'S008',
'S009',
'S010',
'S011',
'S012',
'S013',
'S014',
'S015',
'S016',
'S017',
'S018',
'S019',
'S020',
'S021',
'S022',
'S023',
'S024',
]

ll = []
for ss in subject_inds:
 ll.append(list_all[ss])

print "ll:"
print ll



# Labels/ROIs
#labellist = ['atlleft-lh']#, 'atlright-rh'
tmin, tmax = -0.5, 0.7
event_ids = {'cncrt_wrd': 1, 'abs_wrd': 2}
# frequency bands with variable number of cycles for wavelets
stim_delay = 0.034 # delay in s
"""
frequencies = np.arange(8, 45, 1)
n_cycles = frequencies / float(7)
n_cycles[frequencies<=20] = 2
n_cycles[frequencies<=20] += np.arange(0,1,0.077)
    # n_cycles[np.logical_and((frequencies>10), (frequencies<=20))] = 3
"""

##loop across subjects...
for cnt, meg in enumerate(ll):
    subject=subjects[subject_inds[cnt]]
    print meg
    fname = data_path + meg + 'block_milk_tsss_filt_pnt1_48_raw.fif'
    
    raw = mne.io.Raw(fname, preload=True)
    print "raw loaded"
    include = []  # or stim channels ['STI 014']
    
    ################################################################################
    # jobs and runtime performance
    n_jobs = 1
    
   raw.filter(l_freq=0.1, h_freq=30,  l_trans_bandwidth=0.05,h_trans_bandwidth=0.1, picks=picks, method='fft',filter_length='200s')#l_trans_bandwidth=0.1,
  # 1Hz high pass is often helpful for fitting ICA
    
    picks_meg = mne.pick_types(raw.info, meg=True, eeg=False, eog=False,
                               stim=False, exclude='bads')
    
    ###############################################################################
    # Before applying artifact correction please learn about your actual artifacts
    # by reading :ref:`tut_artifacts_detect`.
    
    ###############################################################################
    # Fit ICA
    # -------
    #
    # ICA parameters:
    
    n_components = 25  # if float, select n_components by explained variance of PCA
    method = 'fastica'  # for comparison with EEGLAB try "extended-infomax" here
    decim = 3  # we need sufficient statistics, not all time points -> saves time
    
    # we will also set state of the random number generator - ICA is a
    # non-deterministic algorithm, but we want to have the same decomposition
    # and the same order of components each time this tutorial is run
    random_state = 23

    ###############################################################################
    # Define the ICA object instance
    ica = ICA(n_components=n_components, method=method, random_state=random_state)
    print(ica)
    
    ###############################################################################
    # we avoid fitting ICA on crazy environmental artifacts that would
    # dominate the variance and decomposition
    reject = dict(mag=5e-12, grad=4000e-13)
    ica.fit(raw, picks=picks_meg, decim=decim, reject=reject)
    print(ica)
    
    ###############################################################################
    # Plot ICA components
    
    ica.plot_components()  # can you spot some potential bad guys?
    
    
    ###############################################################################
    # Component properties
    # --------------------
    #
    # Let's take a closer look at properties of first three independent components.
    
    # first, component 0:
    ica.plot_properties(raw, picks=0)
    
    ###############################################################################
    # we can see that the data were filtered so the spectrum plot is not
    # very informative, let's change that:
    ica.plot_properties(raw, picks=0, psd_args={'fmax': 35.})
    
    ###############################################################################
    # we can also take a look at multiple different components at once:
    ica.plot_properties(raw, picks=[1, 2], psd_args={'fmax': 35.})
    
    ###############################################################################
    # Instead of opening individual figures with component properties, we can
    # also pass an instance of Raw or Epochs in ``inst`` arument to
    # ``ica.plot_components``. This would allow us to open component properties
    # interactively by clicking on individual component topomaps. In the notebook
    # this woks only when running matplotlib in interactive mode (``%matplotlib``).
    
    # uncomment the code below to test the inteactive mode of plot_components:
    # ica.plot_components(picks=range(10), inst=raw)
    
    ###############################################################################
    # Advanced artifact detection
    # ---------------------------
    #
    # Let's use a more efficient way to find artefacts
    
    eog_average = create_eog_epochs(raw, reject=dict(mag=5e-12, grad=4000e-13),
                                    picks=picks_meg).average()
    
    # We simplify things by setting the maximum number of components to reject
    n_max_eog = 1  # here we bet on finding the vertical EOG components
    eog_epochs = create_eog_epochs(raw, reject=reject)  # get single EOG trials
    eog_inds, scores = ica.find_bads_eog(eog_epochs)  # find via correlation
    
    ica.plot_scores(scores, exclude=eog_inds)  # look at r scores of components
    # we can see that only one component is highly correlated and that this
    # component got detected by our correlation analysis (red).
    
    ica.plot_sources(eog_average, exclude=eog_inds)  # look at source time course
    
    ###############################################################################
    # We can take a look at the properties of that component, now using the
    # data epoched with respect to EOG events.
    # We will also use a little bit of smoothing along the trials axis in the
    # epochs image:
    ica.plot_properties(eog_epochs, picks=eog_inds, psd_args={'fmax': 35.},
                        image_args={'sigma': 1.})
    
    ###############################################################################
    # That component is showing a prototypical average vertical EOG time course.
    #
    # Pay attention to the labels, a customized read-out of the
    # :attr:`ica.labels_ <mne.preprocessing.ICA.labels_>`
    print(ica.labels_)
    
    ###############################################################################
    # These labels were used by the plotters and are added automatically
    # by artifact detection functions. You can also manually edit them to annotate
    # components.
    #
    # Now let's see how we would modify our signals if we removed this component
    # from the data
    ica.plot_overlay(eog_average, exclude=eog_inds, show=False)
    # red -> before, black -> after. Yes! We remove quite a lot!
    
    # to definitely register this component as a bad one to be removed
    # there is the ``ica.exclude`` attribute, a simple Python list
    ica.exclude.extend(eog_inds)
    
    # from now on the ICA will reject this component even if no exclude
    # parameter is passed, and this information will be stored to disk
    # on saving
    
    # uncomment this for reading and writing
    # ica.save('my-ica.fif')
    # ica = read_ica('my-ica.fif')
    
    ###############################################################################
    # Exercise: find and remove ECG artifacts using ICA!
    ecg_epochs = create_ecg_epochs(raw, tmin=-.5, tmax=.5)
    ecg_inds, scores = ica.find_bads_ecg(ecg_epochs, method='ctps')
    ica.plot_properties(ecg_epochs, picks=ecg_inds, psd_args={'fmax': 35.})
    
    

