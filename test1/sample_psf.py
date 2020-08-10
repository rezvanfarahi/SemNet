"""
==========================================================
Compute point-spread functions (PSFs) for MNE/dSPM/sLORETA
==========================================================

PSFs are computed for four labels in the MNE sample data set
for linear inverse operators (MNE, dSPM, sLORETA).
PSFs describe the spread of activation from one label
across the cortical surface.
"""

# Authors: Olaf Hauk <olaf.hauk@mrc-cbu.cam.ac.uk>
#          Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
#
# License: BSD (3-clause)



# Author: Martin Luessi <mluessi@nmr.mgh.harvard.edu>
#
# License: BSD (3-clause)

print(__doc__)
print "hi! this script is for spectral functional connectivity. Very nice maps are fruits which worth waiting for I bet!"
import sys

sys.path.append('/imaging/local/software/python_packages/nibabel/1.3.0')
sys.path.append('/imaging/local/software/python_packages/pysurfer/v0.4')
#sys.path.append('/imaging/local/software/python_packages/scikit-learn/v0.14.1')
#sys.path.append('/imaging/local/software/mne_python/latest')

# for qsub
# (add ! here if needed) 
sys.path.append('/imaging/local/software/EPD/latest/x86_64/bin/python')
sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
###

import os
import numpy as np
import mne
from mne import io
from mne.io import Raw
import operator
from mne.minimum_norm import (read_inverse_operator, make_inverse_operator, write_inverse_operator, apply_inverse, apply_inverse_epochs,point_spread_function)
from mne.minimum_norm.inverse import (prepare_inverse_operator)
import sklearn
import scipy.io
from mne import find_events
from mne.connectivity import seed_target_indices, spectral_connectivity
data_path = '/imaging/rf02/TypLexMEG/' # root directory for your MEG data
os.chdir(data_path)
event_path = '/imaging/rf02/TypLexMEG/'    # where event files are
label_path = '/group/erp/data/olaf.hauk/MEG/TypLexM/labels/createdlabels_SL/'
inv_path = '/imaging/rf02/TypLexMEG/'
subjects_dir = data_path 
print sys.argv
subject_inds = [1]
for ss in sys.argv[1:]:
   subject_inds.append( int( ss ) )

#subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
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
'meg11_0147/110603/'
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
ii=-1
for meg in ll:
	ii=ii+1
	print ii
	subjects_dir = data_path 
	fname_fwd = data_path + meg + 'forward_5-3L-EMEG-fwd.fif'
	fname_inv = inv_path + meg + 'InverseOperator_fft_1_48_clean_ica_EMEG-inv.fif'
	fname_label = [label_path + subjects[subject_inds[ii]] + '/' + 'atlleft-lh.label',
		       label_path + subjects[subject_inds[ii]] + '/' + 'atlright-rh.label',
		       ]


	# read forward solution (sources in surface-based coordinates)
	forward = mne.read_forward_solution(fname_fwd, force_fixed=False,
		                            surf_ori=True)

	# read inverse operators
	inverse_operator_eegmeg = read_inverse_operator(fname_inv)
	#inverse_operator_meg = read_inverse_operator(fname_inv_meg)

	# read label(s)
	labels = [mne.read_label(ss) for ss in fname_label]

	# regularisation parameter
	snr = 3.0
	lambda2 = 1.0 / snr ** 2
	method = 'MNE'  # can be 'MNE' or 'sLORETA'
	mode = 'svd'
	n_svd_comp = 1

	stc_psf_eegmeg, _ = point_spread_function(inverse_operator_eegmeg, forward, method=method,  labels=labels, lambda2=lambda2, pick_ori='normal', mode=mode, n_svd_comp=n_svd_comp)

	#stc_psf_meg, _ = point_spread_function(inverse_operator_meg, forward, method=method,labels=labels, lambda2=lambda2, pick_ori='normal', mode=mode,n_svd_comp=n_svd_comp)

	# save for viewing in mne_analyze in order of labels in 'labels'
	# last sample is average across PSFs
	data_out = data_path + meg +  'psf_eegmeg'
	subject_from = subjects[subject_inds[ii]]
    	subject_to = 'avgsubj'
    	vertices_to = [np.arange(10242), np.arange(10242)]
	stc_to = mne.morph_data(subject_from, subject_to, stc_psf_eegmeg, subjects_dir=subjects_dir, n_jobs=1, grade=vertices_to)
	stc_to.save(data_out)
	# stc_psf_meg.save('psf_meg')
	"""
	from mayavi import mlab
	fmin = 0.
	time_label = "EEGMEG %d"
	fmax = stc_psf_eegmeg.data[:, 0].max()
	fmid = fmax / 2.
	brain_eegmeg = stc_psf_eegmeg.plot(surface='inflated', hemi='rh',
		                           subjects_dir=subjects_dir,
		                           time_label=time_label, fmin=fmin,
		                           fmid=fmid, fmax=fmax,
		                           figure=mlab.figure(size=(500, 500)))

	time_label = "MEG %d"
	fmax = stc_psf_meg.data[:, 0].max()
	fmid = fmax / 2.
	brain_meg = stc_psf_meg.plot(surface='inflated', hemi='rh',
		                     subjects_dir=subjects_dir,
		                     time_label=time_label, fmin=fmin,
		                     fmid=fmid, fmax=fmax,
		                     figure=mlab.figure(size=(500, 500)))

	# The PSF is centred around the right auditory cortex label,
	# but clearly extends beyond it.
	# It also contains "sidelobes" or "ghost sources"
	# in middle/superior temporal lobe.
	# For the Aud-RH example, MEG and EEGMEG do not seem to differ a lot,
	# but the addition of EEG still decreases point-spread to distant areas
	# (e.g. to ATL and IFG).
	# The chosen labels are quite far apart from each other, so their PSFs
	# do not overlap (check in mne_analyze)
	"""
