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
sys.path.append('/imaging/local/software/anaconda/latest/x86_64/bin/python')
sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
###

import matplotlib.pyplot as plt
import numpy as np
import mne
from mne import io
from mne.io import Raw
import pylab as pl
from mne import (read_evokeds, equalize_channels,read_proj, read_selection)
from mne.preprocessing.ica import ICA
import scipy.io as sio
import operator

from mne.minimum_norm import (read_inverse_operator, make_inverse_operator, write_inverse_operator, apply_inverse, read_inverse_operator,apply_inverse_epochs, apply_inverse, source_induced_power,source_band_induced_power)
from mne.minimum_norm.inverse import (prepare_inverse_operator, _assemble_kernel)

from surfer import Brain
from surfer.io import read_stc
import logging
import os
import sklearn
import scipy.io
from mne import filter
from mne import find_events
from mne.epochs import combine_event_ids

from mne.layouts import read_layout
from mne.time_frequency import compute_raw_psd
from mne.connectivity import seed_target_indices, spectral_connectivity

###############################################################################
data_path = '/imaging/rf02/TypLexMEG/' # root directory for your MEG data
out_path = '/imaging/rf02/TypLexMEG/icaanalysis_results/stc/GrandAverage/connectivity/'
subjects_dir = '/imaging/rf02/TypLexMEG/'    # where your MRI subdirectories are
# where event-files are
event_path = '/imaging/rf02/TypLexMEG/'    # where event files are

label_path = '/home/rf02/rezvan/TypLexMEG/createdlabels_SL/'


# get indices for subjects to be processed from command line input
# 
print sys.argv
subject_inds = []
for ss in sys.argv[1:]:
   subject_inds.append( int( ss ) )

subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
print "subject_inds:"
print subject_inds
print "No rejection"

list_all = ['meg10_0378/101209/', 
'meg10_0390/101214/',
'meg11_0026/110223/', 
'meg11_0050/110307/', 
'meg11_0052/110307/', 
'meg11_0069/110315/', 
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
'meg11_0147/110603/', 


]

# subjects names used for MRI data
subjects=['SS1_Meg10_0378',
'SS2_Meg10_0390',
'SS3_Meg11_0026',
'SS4_Meg11_0050',
'SS5_Meg11_0052',
'SS6_Meg11_0069',
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
#labellist = ['atlleft-lh']#, 'atlleftt-rh'
tmin, tmax = -0.5, 0.7
reject = dict(eeg=120e-6,grad=200e-12, mag=4e-12)
event_ids = {'cncrt_wrd': 1, 'abs_wrd': 2}
# frequency bands with variable number of cycles for wavelets
stim_delay = 0.034 # delay in s
frequencies = np.arange(8, 45, 1)
n_cycles = frequencies / float(7)
n_cycles[frequencies<=20] = 2
n_cycles[frequencies<=20] += np.arange(0,1,0.077)
    # n_cycles[np.logical_and((frequencies>10), (frequencies<=20))] = 3
labellist = ['atlleft-lh', 'atlright-rh', 'medialtempright-rh','medialtempleft-lh']

for label_name in labellist:
	"""
	fname0 = data_path + ll[0] + 'Morphed_SemLoc_icaclean_Abstract_PhaseLocking_ThetaAlphaBetaGamma_150_350_' + label_name[0:-3]
	stc0 = mne.read_source_estimate(fname0)
	

	fname1 = data_path + ll[1] + 'Morphed_SemLoc_icaclean_Abstract_PhaseLocking_ThetaAlphaBetaGamma_150_350_' + label_name[0:-3]
	stc1 = mne.read_source_estimate(fname1)

	fname2 = data_path + ll[2] + 'Morphed_SemLoc_icaclean_Abstract_PhaseLocking_ThetaAlphaBetaGamma_150_350_' + label_name[0:-3]
	stc2 = mne.read_source_estimate(fname2)

	fname3 = data_path + ll[3] + 'Morphed_SemLoc_icaclean_Abstract_PhaseLocking_ThetaAlphaBetaGamma_150_350_' + label_name[0:-3]
	stc3 = mne.read_source_estimate(fname3)

	fname4 = data_path + ll[4] + 'Morphed_SemLoc_icaclean_Abstract_PhaseLocking_ThetaAlphaBetaGamma_150_350_' + label_name[0:-3]
	stc4 = mne.read_source_estimate(fname4)

	fname5 = data_path + ll[5] + 'Morphed_SemLoc_icaclean_Abstract_PhaseLocking_ThetaAlphaBetaGamma_150_350_' + label_name[0:-3]
	stc5 = mne.read_source_estimate(fname5)

	fname6 = data_path + ll[6] + 'Morphed_SemLoc_icaclean_Abstract_PhaseLocking_ThetaAlphaBetaGamma_150_350_' + label_name[0:-3]
	stc6 = mne.read_source_estimate(fname6)

	fname7 = data_path + ll[7] + 'Morphed_SemLoc_icaclean_Abstract_PhaseLocking_ThetaAlphaBetaGamma_150_350_' + label_name[0:-3]
	stc7 = mne.read_source_estimate(fname7)

	fname8 = data_path + ll[8] + 'Morphed_SemLoc_icaclean_Abstract_PhaseLocking_ThetaAlphaBetaGamma_150_350_' + label_name[0:-3]
	stc8 = mne.read_source_estimate(fname8)

	fname9 = data_path + ll[9] + 'Morphed_SemLoc_icaclean_Abstract_PhaseLocking_ThetaAlphaBetaGamma_150_350_' + label_name[0:-3]
	stc9 = mne.read_source_estimate(fname9)

	fname10 = data_path + ll[10] + 'Morphed_SemLoc_icaclean_Abstract_PhaseLocking_ThetaAlphaBetaGamma_150_350_' + label_name[0:-3]
	stc10 = mne.read_source_estimate(fname10)

	fname11 = data_path + ll[11] + 'Morphed_SemLoc_icaclean_Abstract_PhaseLocking_ThetaAlphaBetaGamma_150_350_' + label_name[0:-3]
	stc11 = mne.read_source_estimate(fname11)

	fname12 = data_path + ll[12] + 'Morphed_SemLoc_icaclean_Abstract_PhaseLocking_ThetaAlphaBetaGamma_150_350_' + label_name[0:-3]
	stc12 = mne.read_source_estimate(fname12)

	fname13 = data_path + ll[13] + 'Morphed_SemLoc_icaclean_Abstract_PhaseLocking_ThetaAlphaBetaGamma_150_350_' + label_name[0:-3]
	stc13 = mne.read_source_estimate(fname13)

	fname14 = data_path + ll[14] + 'Morphed_SemLoc_icaclean_Abstract_PhaseLocking_ThetaAlphaBetaGamma_150_350_' + label_name[0:-3]
	stc14 = mne.read_source_estimate(fname14)

	fname15 = data_path + ll[15] + 'Morphed_SemLoc_icaclean_Abstract_PhaseLocking_ThetaAlphaBetaGamma_150_350_' + label_name[0:-3]
	stc15 = mne.read_source_estimate(fname15)

	fname16 = data_path + ll[16] + 'Morphed_SemLoc_icaclean_Abstract_PhaseLocking_ThetaAlphaBetaGamma_150_350_' + label_name[0:-3]
	stc16 = mne.read_source_estimate(fname16)

	stc_grand=np.mean([stc0,stc1,stc2,stc3,stc4,stc5,stc6,stc7,stc8,stc9,stc10,stc11,stc12,stc13,stc14,stc15,stc16])#

	out_file=out_path + 'GrandAverage_SemLoc_icaclean_Abstract_PhaseLocking_ThetaAlphaBetaGamma_150_350_' + label_name[0:-3]
	stc_grand.save(out_file)
	
	"""
	out_file1=out_path + 'GrandAverage_SemLoc_icaclean_Concrete_PhaseLocking_ThetaAlphaBetaGamma_150_350_' + label_name[0:-3]
	stcs1 = mne.read_source_estimate(out_file1)
	out_file2=out_path + 'GrandAverage_SemLoc_icaclean_Abstract_PhaseLocking_ThetaAlphaBetaGamma_150_350_' + label_name[0:-3]
	stcs2 = mne.read_source_estimate(out_file2)
	stc_subtract=np.subtract(stcs1,stcs2)
	out_file=out_path + 'GrandAverage_SemLoc_icaclean_Subtract_PhaseLocking_ThetaAlphaBetaGamma_150_350_' + label_name[0:-3]
	stc_subtract.save(out_file)
	
	
"""
plot(subject='average', surface='inflated', hemi='splited', colormap='hot', time_label='time=%0.2f ms', smoothing_steps=5, fmin=0, fmid=.25, fmax=.5, transparent=True, alpha=1.0, time_viewer=False, config_opts={}, subjects_dir=subjects_dir, figure=None, views='lat', colorbar=True)
brain.save_image('DICS_source_power_freq_%d.png' % csd.frequencies[0])
"""
