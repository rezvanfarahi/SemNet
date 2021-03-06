"""
=========================================================
TF representation of SemLoc for frequency bands in source labels
=========================================================

"""
# Authors: Alexandre Gramfort <gramfort@nmr.mgh.harvard.edu>
#
# License: BSD (3-clause)


print(__doc__)

import sys
sys.path.append('/imaging/local/software/python_packages/nibabel/1.3.0')
sys.path.append('/imaging/local/software/python_packages/pysurfer/v0.4')
sys.path.append('/imaging/local/software/python_packages/scikit-learn/v0.14.1')
sys.path.append('/imaging/local/software/mne_python/latest')

# for qsub
# (add ! here if needed) 
#sys.path.append('/imaging/local/software/anaconda/latest/x86_64/bin/python')
sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
###
import os
print "hi"
#import matplotlib.pyplot as plt
import numpy as np
print "hi1"
import mne
import scipy
print "hi2"
from mne import io
print "hi3"
from mne.io import Raw
print "hi4"
#import pylab as pl
#from mne import (fiff, read_evokeds, equalize_channels,read_proj, read_selection)
#from mne.preprocessing.ica import ICA
import scipy.io as sio
print "hi5"
import operator
print "hi6"

from mne.minimum_norm import (read_inverse_operator, make_inverse_operator, write_inverse_operator, apply_inverse, apply_inverse_epochs, source_induced_power,source_band_induced_power)
print "hi7"
from mne.minimum_norm.inverse import (prepare_inverse_operator)
print "hi8"

#from surfer import Brain
#from surfer.io import read_stc
#import logging

import sklearn
print "hi9"
import scipy.io
print "hi10"
#from mne import filter
from mne import find_events
print "hi11"
#from mne.epochs import combine_event_ids
#from mne.layouts import read_layout
from mne.time_frequency import compute_raw_psd, induced_power
print "hi12"
#from mne.connectivity import seed_target_indices, spectral_connectivity

###############################################################################
data_path = '/imaging/rf02/TypLexMEG/' # root directory for your MEG data
#os.chdir(data_path)
subjects_dir = '/imaging/rf02/TypLexMEG/'    # where your MRI subdirectories are
# where event-files are
event_path = '/imaging/rf02/TypLexMEG/'    # where event files are

label_path = '/imaging/rf02/TypLexMEG/createdlabels_SL/'

inv_path = '/imaging/rf02/TypLexMEG/'
inv_fname = 'InverseOperator_EMEG-inv.fif'

# get indices for subjects to be processed from command line input
# 
print sys.argv
subject_inds = []
for ss in sys.argv[1:]:
   subject_inds.append( int( ss ) )

subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
print "subject_inds:"
print subject_inds
print "No rejection"

list_all = ['meg10_0378/101209/', 
'meg10_0390/101214/',
'meg11_0026/110223/', 
'meg11_0050/110307/', 
'meg11_0052/110307/', 
#'meg11_0069/110315/', 
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

# subjects names used for MRI data
subjects=['SS1_Meg10_0378',
'SS2_Meg10_0390',
'SS3_Meg11_0026',
'SS4_Meg11_0050',
'SS5_Meg11_0052',
#'SS6_Meg11_0069',
'SS9_Meg11_0086',
'SS10_Meg11_0091',
'SS12_Meg11_0096',
'SS13_Meg11_0101',
'SS14_Meg11_0102',
'SS15_Meg11_0112',
'SS16_Meg11_0104',
'SS18_Meg11_0118',
'SS19_Meg11_0131',
'SS20_Meg11_0144',
'SS21_Meg11_0147'
]

ll = []
for ss in subject_inds:
 ll.append(list_all[ss])

print "ll:"
print ll



# Labels/ROIs
labellist = ['atlleft-lh'] #['atlleft-lh',
tmin, tmax = -0.5, 0.7
reject = dict(eeg=120e-6, eog=150e-6, grad=200e-12, mag=4e-12)
event_ids = {'cncrt_wrd': 1, 'abs_wrd': 2}
# frequency bands with variable number of cycles for wavelets
stim_delay = 0.034 # delay in s
bands=['gamma']#,'alpha','beta', 'gamma']
"""
frequencies = np.arange(8, 45, 1)
n_cycles = frequencies / float(7)
n_cycles[frequencies<=20] = 2
n_cycles[frequencies<=20] += np.arange(0,1,0.077)
    # n_cycles[np.logical_and((frequencies>10), (frequencies<=20))] = 3
"""

bb1=np.zeros((16,1))
bb2=np.zeros((16,1))
bb3=np.zeros((16,1))
bb4=np.zeros((16,1))
bb5=np.zeros((16,1))
bb6=np.zeros((16,1))
bb7=np.zeros((16,1))
##loop across subjects...
for ii, meg in enumerate(ll):
	kk=-1
	print ii
	for b in bands:
		print b
		fname1 = data_path + meg + 'Morphed_SemLoc_Concrete_Power_m500_700' + b  
		stc1 = mne.read_source_estimate(fname1)
		if b=='theta':
 		   fname2 = data_path + meg + 'Morphed_SemLoc_Abstract_Power_m500_700_' + b  
		else:
		   fname2 = data_path + meg + 'Morphed_SemLoc_Abstract_Power_m500_700' + b 
		stc2 = mne.read_source_estimate(fname2)
		stc=np.subtract(stc1,stc2)
		for label_name in labellist:
			kk=kk+1
			print kk
			fname_label = label_path + label_name + '.label'
			label1 = mne.read_label(fname_label)
			stc_label = stc.in_label(label1) 
			bb1[ii,kk]=stc_label.data[:,590:610].mean()
			bb2[ii,kk]=stc_label.data[:,660:700].mean()
			bb3[ii,kk]=stc_label.data[:,720:760].mean()
			bb4[ii,kk]=stc_label.data[:,820:860].mean()
			bb5[ii,kk]=stc_label.data[:,900:1100].mean()
			bb6[ii,kk]=stc_label.data[:,650:850].mean()
			bb7[ii,kk]=stc_label.data[:,750:950].mean()

#tvals1, cvals1, pvals1, h0val1 = mne.stats.permutation_cluster_1samp_test(bb1, threshold=None, n_permutations=1024, tail=0, connectivity=None, verbose=None, n_jobs=1)
t1, p1, h01 = mne.stats.permutation_t_test(bb1, n_permutations=500, tail=0, n_jobs=1, verbose=None)
t2, p2, h02 = mne.stats.permutation_t_test(bb2, n_permutations=500, tail=0, n_jobs=1, verbose=None)
t3, p3, h03 = mne.stats.permutation_t_test(bb3, n_permutations=500, tail=0, n_jobs=1, verbose=None)
t4, p4, h04 = mne.stats.permutation_t_test(bb4, n_permutations=500, tail=0, n_jobs=1, verbose=None)
t5, p5, h05 = mne.stats.permutation_t_test(bb5, n_permutations=500, tail=0, n_jobs=1, verbose=None)
t6, p6, h06 = mne.stats.permutation_t_test(bb6, n_permutations=500, tail=0, n_jobs=1, verbose=None)
t7, p7, h07 = mne.stats.permutation_t_test(bb7, n_permutations=500, tail=0, n_jobs=1, verbose=None)
print t1,p1; print t2,p2; print t3,p3; print t4,p4; print t5,p5; print t6,p6; print t7,p7 
#t1 = mne.stats.ttest_1samp_no_p(bb1, sigma=1e-3, method='relative')
#t2 = mne.stats.ttest_1samp_no_p(bb2, sigma=1e-3, method='relative')
#t3 = mne.stats.ttest_1samp_no_p(bb3, sigma=1e-3, method='relative')
#t4 = mne.stats.ttest_1samp_no_p(bb4, sigma=1e-3, method='relative')
#t5 = mne.stats.ttest_1samp_no_p(bb5, sigma=1e-3, method='relative')
#t6 = mne.stats.ttest_1samp_no_p(bb6, sigma=1e-3, method='relative')
#t7 = mne.stats.ttest_1samp_no_p(bb7, sigma=1e-3, method='relative')
p_threshold = 0.05
n_subjects=16
t_threshold = -scipy.stats.distributions.t.ppf(p_threshold / 2., n_subjects - 1)

			
			#fname_out = data_path + meg + 'Morphed_SemLoc_Abstract_Power_m500_700_' + b + label_name
			#stc_label.save(fname_out) 

    
