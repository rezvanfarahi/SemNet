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
import os.path as op
import matplotlib.pyplot as plt
import numpy as np
import mne
from mne import io
from mne.io import Raw
import pylab as pl
from mne import (stats, read_evokeds, equalize_channels,read_proj, read_selection)
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
from mne import (io, spatial_tris_connectivity, compute_morph_matrix,
                 grade_to_tris)
from mne.epochs import equalize_epoch_counts
from mne.stats import (spatio_temporal_cluster_1samp_test,
                       summarize_clusters_stc)
from mne.minimum_norm import apply_inverse, read_inverse_operator
from mne.datasets import sample
from mne.fiff import Evoked
from mne.layouts import read_layout
from mne.time_frequency import compute_raw_psd, induced_power
from mne.connectivity import seed_target_indices, spectral_connectivity

###############################################################################
data_path = '/imaging/rf02/TypLexMEG/' # root directory for your MEG data
os.chdir(data_path)
subjects_dir = '/home/rf02/rezvan/TypLexMEG/'    # where your MRI subdirectories are
# where event-files are
event_path = '/imaging/rf02/TypLexMEG/'    # where event files are

label_path = '/home/rf02/rezvan/TypLexMEG/createdlabels_SL/'

inv_path = '/imaging/rf02/TypLexMEG/'
inv_fname = 'InverseOperator_EMEG-inv.fif'

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
tmin, tmax = -0.3, 0.6
reject = dict(eeg=120e-6, eog=150e-6, grad=200e-12, mag=4e-12)
event_ids = {'cncrt_wrd': 1, 'abs_wrd': 2}
# frequency bands with variable number of cycles for wavelets
stim_delay = 0.034 # delay in s
frequencies = np.arange(8, 45, 1)
n_cycles = frequencies / float(7)
n_cycles[frequencies<=20] = 2
n_cycles[frequencies<=20] += np.arange(0,1,0.077)
    # n_cycles[np.logical_and((frequencies>10), (frequencies<=20))] = 3
labellist = ['atlleft-lh', 'atlright-rh','medialtempright-rh','medialtempleft-lh']


for label_name in labellist:
	ii=-1
	#X=np.zeros((20484,2,17))
	X_theta=np.zeros((20484,16))
	X_alpha=np.zeros((20484,16))
	X_beta=np.zeros((20484,16))
	X_gamma=np.zeros((20484,16))
	for meg in ll:
		ii=ii+1;
		print ii
	
### Concrete
		fname_cncrt = data_path + meg + 'Morphed_SemLoc_Concrete_Coherence_AlphaBeta_250_450_' + label_name[0:-3]
		stc_cncrt = mne.read_source_estimate(fname_cncrt)
		fname_abs = data_path + meg + 'Morphed_SemLoc_Abstract_Coherence_AlphaBeta_250_450_' + label_name[0:-3]
		stc_abs = mne.read_source_estimate(fname_abs)	
		stc_subtract=np.subtract(stc_cncrt,stc_abs)
		#X[:,:,ii]=stc_subtract.data
		X_theta[:,ii]=stc_subtract.data[:,0]
		X_alpha[:,ii]=stc_subtract.data[:,1]
		X_beta[:,ii]=stc_subtract.data[:,2]
		X_gamma[:,ii]=stc_subtract.data[:,3]	
		
	#X = np.transpose(X, [2, 1, 0])	
	connectivity = spatial_tris_connectivity(grade_to_tris(5))
	threshold = 0.5
	fsave_vertices = [np.arange(10242), np.arange(10242)]
	n_permutations=1024
	X_theta=np.transpose(X_theta,[1,0])
	X_alpha=np.transpose(X_alpha,[1,0])
	X_beta=np.transpose(X_beta,[1,0])
	X_gamma=np.transpose(X_gamma,[1,0])
	
	t_theta= mne.stats.ttest_1samp_no_p(X_theta, sigma=0, method='relative')
	print "per1 done"
	t_alpha = mne.stats.ttest_1samp_no_p(X_alpha, sigma=0, method='relative')
	print "per2 done"
	t_beta = mne.stats.ttest_1samp_no_p(X_beta, sigma=0, method='relative')
	print "per3 done"
	t_gamma = mne.stats.ttest_1samp_no_p(X_gamma, sigma=0, method='relative')
	print "per4 done"

	#t_map=np.concatenate((t_theta,t_alpha,t_beta,t_gamma)).reshape((20484,4))
	t_map=np.zeros((20484,4))
	t_map[:,0]=t_theta
	t_map[:,1]=t_alpha
	t_map[:,2]=t_beta
	t_map[:,3]=t_gamma
	#p_map=np.concatenate((p_theta,p_alpha,p_beta,p_gamma)).reshape((20484,4))
	t_stc = mne.SourceEstimate(t_map, vertices=fsave_vertices, tmin=1e-3*2, tstep=1e-3*2, subject='average', verbose=None)
	print "se1 done"
	#p_stc = mne.SourceEstimate(p_map, vertices=fsave_vertices, tmin=1e-3*2, tstep=1e-3*2, subject='average', verbose=None)
	print "se1 done"
	#t_beta_stc = mne.SourceEstimate(t_beta, tmin=0, tstep=0, vertices=fsave_vertices)
	#p_beta_stc = mne.SourceEstimate(p_beta, tmin=0, tstep=0, vertices=fsave_vertices)
	"""
	X1=np.zeros((20484,4))
	X1[:,0]=X_theta[0,:]
	X1[:,1]=X_alpha[0,:]
	X1[:,2]=X_beta[0,:]
	X1[:,3]=X_gamma[0,:]
	Xstc=t_stc = mne.SourceEstimate(X1, vertices=fsave_vertices, tmin=1e-3*2, tstep=1e-3*2, subject='average', verbose=None)
	out_filet=data_path + 'UVttest_Test' + label_name[0:-3]
	t_stc.save(out_filet)
	"""
	
	"""
	print('Computing connectivity.')
	connectivity = spatial_tris_connectivity(grade_to_tris(5))

	#    Note that X needs to be a multi-dimensional array of shape
	#    samples (subjects) x time x space, so we permute dimensions
	X = np.transpose(X, [2, 1, 0])

	#    Now let's actually do the clustering. This can take a long time...
	#    Here we set the threshold quite high to reduce computation.
	p_threshold = 0.001
	n_subjects=17
	#t_threshold = -stats.distributions.t.ppf(p_threshold / 2., n_subjects - 1)
	print('Clustering.')
	T_obs, clusters, cluster_p_values, H0 = clu = spatio_temporal_cluster_1samp_test(X, connectivity=connectivity, n_jobs=2, threshold=None)
	#    Now select the clusters that are sig. at p < 0.05 (note that this value
	#    is multiple-comparisons corrected).
	#good_cluster_inds = np.where(cluster_p_values < 0.05)[0]
	tstep=stc_abs.tstep
	fsave_vertices = [np.arange(10242), np.arange(10242)]
	stc_all_cluster_vis = summarize_clusters_stc(clu, tstep=tstep,
                                             vertno=fsave_vertices,
                                             subject='average')
	"""
	out_file1=data_path + 'UVttest_tmap_SemLoc_ConcreteAbstract_Coherence_ThetaAlphaBetaGamma_250_450_' + label_name[0:-3]
	t_stc.save(out_file1)
	#out_file2=data_path + 'Permutation_pmap_SemLoc_ConcreteAbstract_Coherence_ThetaAlphaBetaGamma_250_450_' + label_name[0:-3]
	#p_stc.save(out_file2)
	#out_file3=data_path + 'Permutation_tmap_SemLoc_ConcreteAbstract_Coherence_Beta_250_450_' + label_name[0:-3]
	#t_beta_stc.save(out_file3)
	#out_file4=data_path + 'Permutation_pmap_SemLoc_ConcreteAbstract_Coherence_Beta_250_450_' + label_name[0:-3]
	#p_beta_stc.save(out_file4)

	#brain = t_stc.plot('average', 'inflated', 'rh', fmin=0., fmid=2., fmax=4., time_label='tmap',subjects_dir=subjects_dir, time_viewer=False, colorbar=True)
	#fig_out = data_path + label_name[0:-3] + '_TEST_coherence.jpg'
	#brain.show_view('lateral')
	#brain.save_image(fig_out)


"""
		fname_cncrt1 = data_path + ll[1] + 'Morphed_SemLoc_Concrete_Coherence_AlphaBeta_250_450_' + label_name[0:-3]
		stc_cncrt1 = mne.read_source_estimate(fname_cncrt1)

		fname_cncrt2 = data_path + ll[2] + 'Morphed_SemLoc_Concrete_Coherence_AlphaBeta_250_450_' + label_name[0:-3]
		stc_cncrt2 = mne.read_source_estimate(fname_cncrt2)

		fname_cncrt3 = data_path + ll[3] + 'Morphed_SemLoc_Concrete_Coherence_AlphaBeta_250_450_' + label_name[0:-3]
		stc_cncrt3 = mne.read_source_estimate(fname_cncrt3)

		fname_cncrt4 = data_path + ll[4] + 'Morphed_SemLoc_Concrete_Coherence_AlphaBeta_250_450_' + label_name[0:-3]
		stc_cncrt4 = mne.read_source_estimate(fname_cncrt4)

		fname_cncrt5 = data_path + ll[5] + 'Morphed_SemLoc_Concrete_Coherence_AlphaBeta_250_450_' + label_name[0:-3]
		stc_cncrt5 = mne.read_source_estimate(fname_cncrt5)

		fname_cncrt6 = data_path + ll[6] + 'Morphed_SemLoc_Concrete_Coherence_AlphaBeta_250_450_' + label_name[0:-3]
		stc_cncrt6 = mne.read_source_estimate(fname_cncrt6)

		fname_cncrt7 = data_path + ll[7] + 'Morphed_SemLoc_Concrete_Coherence_AlphaBeta_250_450_' + label_name[0:-3]
		stc_cncrt7 = mne.read_source_estimate(fname_cncrt7)

		fname_cncrt8 = data_path + ll[8] + 'Morphed_SemLoc_Concrete_Coherence_AlphaBeta_250_450_' + label_name[0:-3]
		stc_cncrt8 = mne.read_source_estimate(fname_cncrt8)

		fname_cncrt9 = data_path + ll[9] + 'Morphed_SemLoc_Concrete_Coherence_AlphaBeta_250_450_' + label_name[0:-3]
		stc_cncrt9 = mne.read_source_estimate(fname_cncrt9)

		fname_cncrt10 = data_path + ll[10] + 'Morphed_SemLoc_Concrete_Coherence_AlphaBeta_250_450_' + label_name[0:-3]
		stc_cncrt10 = mne.read_source_estimate(fname_cncrt10)

		fname_cncrt11 = data_path + ll[11] + 'Morphed_SemLoc_Concrete_Coherence_AlphaBeta_250_450_' + label_name[0:-3]
		stc_cncrt11 = mne.read_source_estimate(fname_cncrt11)

		fname_cncrt12 = data_path + ll[12] + 'Morphed_SemLoc_Concrete_Coherence_AlphaBeta_250_450_' + label_name[0:-3]
		stc_cncrt12 = mne.read_source_estimate(fname_cncrt12)

		fname_cncrt13 = data_path + ll[13] + 'Morphed_SemLoc_Concrete_Coherence_AlphaBeta_250_450_' + label_name[0:-3]
		stc_cncrt13 = mne.read_source_estimate(fname_cncrt13)

		fname_cncrt14 = data_path + ll[14] + 'Morphed_SemLoc_Concrete_Coherence_AlphaBeta_250_450_' + label_name[0:-3]
		stc_cncrt14 = mne.read_source_estimate(fname_cncrt14)

		fname_cncrt15 = data_path + ll[15] + 'Morphed_SemLoc_Concrete_Coherence_AlphaBeta_250_450_' + label_name[0:-3]
		stc_cncrt15 = mne.read_source_estimate(fname_cncrt15)

		fname_cncrt16 = data_path + ll[16] + 'Morphed_SemLoc_Concrete_Coherence_AlphaBeta_250_450_' + label_name[0:-3]
		stc_cncrt16 = mne.read_source_estimate(fname_cncrt16)

	### Abstract
		
		fname_abs1 = data_path + ll[1] + 'Morphed_SemLoc_Abstract_Coherence_AlphaBeta_250_450_' + label_name[0:-3]
		stc_abs1 = mne.read_source_estimate(fname_abs1)

		fname_abs2 = data_path + ll[2] + 'Morphed_SemLoc_Abstract_Coherence_AlphaBeta_250_450_' + label_name[0:-3]
		stc_abs2 = mne.read_source_estimate(fname_abs2)

		fname_abs3 = data_path + ll[3] + 'Morphed_SemLoc_Abstract_Coherence_AlphaBeta_250_450_' + label_name[0:-3]
		stc_abs3 = mne.read_source_estimate(fname_abs3)

		fname_abs4 = data_path + ll[4] + 'Morphed_SemLoc_Abstract_Coherence_AlphaBeta_250_450_' + label_name[0:-3]
		stc_abs4 = mne.read_source_estimate(fname_abs4)

		fname_abs5 = data_path + ll[5] + 'Morphed_SemLoc_Abstract_Coherence_AlphaBeta_250_450_' + label_name[0:-3]
		stc_abs5 = mne.read_source_estimate(fname_abs5)

		fname_abs6 = data_path + ll[6] + 'Morphed_SemLoc_Abstract_Coherence_AlphaBeta_250_450_' + label_name[0:-3]
		stc_abs6 = mne.read_source_estimate(fname_abs6)

		fname_abs7 = data_path + ll[7] + 'Morphed_SemLoc_Abstract_Coherence_AlphaBeta_250_450_' + label_name[0:-3]
		stc_abs7 = mne.read_source_estimate(fname_abs7)

		fname_abs8 = data_path + ll[8] + 'Morphed_SemLoc_Abstract_Coherence_AlphaBeta_250_450_' + label_name[0:-3]
		stc_abs8 = mne.read_source_estimate(fname_abs8)

		fname_abs9 = data_path + ll[9] + 'Morphed_SemLoc_Abstract_Coherence_AlphaBeta_250_450_' + label_name[0:-3]
		stc_abs9 = mne.read_source_estimate(fname_abs9)

		fname_abs10 = data_path + ll[10] + 'Morphed_SemLoc_Abstract_Coherence_AlphaBeta_250_450_' + label_name[0:-3]
		stc_abs10 = mne.read_source_estimate(fname_abs10)

		fname_abs11 = data_path + ll[11] + 'Morphed_SemLoc_Abstract_Coherence_AlphaBeta_250_450_' + label_name[0:-3]
		stc_abs11 = mne.read_source_estimate(fname_abs11)

		fname_abs12 = data_path + ll[12] + 'Morphed_SemLoc_Abstract_Coherence_AlphaBeta_250_450_' + label_name[0:-3]
		stc_abs12 = mne.read_source_estimate(fname_abs12)

		fname_abs13 = data_path + ll[13] + 'Morphed_SemLoc_Abstract_Coherence_AlphaBeta_250_450_' + label_name[0:-3]
		stc_abs13 = mne.read_source_estimate(fname_abs13)

		fname_abs14 = data_path + ll[14] + 'Morphed_SemLoc_Abstract_Coherence_AlphaBeta_250_450_' + label_name[0:-3]
		stc_abs14 = mne.read_source_estimate(fname_abs14)

		fname_abs15 = data_path + ll[15] + 'Morphed_SemLoc_Abstract_Coherence_AlphaBeta_250_450_' + label_name[0:-3]
		stc_abs15 = mne.read_source_estimate(fname_abs15)

		fname_abs16 = data_path + ll[16] + 'Morphed_SemLoc_Abstract_Coherence_AlphaBeta_250_450_' + label_name[0:-3]
		stc_abs16 = mne.read_source_estimate(fname_abs16)


		###contrast
		stc_subtract0=np.subtract(stc_cncrt0,stc_abs0)
		stc_subtract1=np.subtract(stc_cncrt1,stc_abs1)
		stc_subtract2=np.subtract(stc_cncrt2,stc_abs2)
		stc_subtract3=np.subtract(stc_cncrt3,stc_abs3)
		stc_subtract4=np.subtract(stc_cncrt4,stc_abs4)
		stc_subtract5=np.subtract(stc_cncrt5,stc_abs5)
		stc_subtract6=np.subtract(stc_cncrt6,stc_abs6)
		stc_subtract7=np.subtract(stc_cncrt7,stc_abs7)
		stc_subtract8=np.subtract(stc_cncrt8,stc_abs8)
		stc_subtract9=np.subtract(stc_cncrt9,stc_abs9)
		stc_subtract10=np.subtract(stc_cncrt10,stc_abs10)
		stc_subtract11=np.subtract(stc_cncrt11,stc_abs11)
		stc_subtract12=np.subtract(stc_cncrt12,stc_abs12)
		stc_subtract13=np.subtract(stc_cncrt13,stc_abs13)
		stc_subtract14=np.subtract(stc_cncrt14,stc_abs14)
		stc_subtract15=np.subtract(stc_cncrt15,stc_abs15)
		stc_subtract16=np.subtract(stc_cncrt16,stc_abs16)
		stc_shape=stc_abs0.data.shape
		X_alpha=np.zeros((stc_shape[0],17))
		for ii in np.arange(17)
		X_alpha[]
		X_beta=np.zeros((stc_shape[0],17))
		[t_alpha,p_alpha,H0_alpha]=permutation_t_test(X_alpha, n_permutations=10000, tail=0, n_jobs=1, verbose=None)

"""


"""
		out_file1=data_path + 'GrandAverage_SemLoc_Concrete_Coherence_AlphaBeta_250_450_' + label_name[0:-3]
		stcs1 = mne.read_source_estimate(out_file1)
		out_file2=data_path + 'GrandAverage_SemLoc_Abstract_Coherence_AlphaBeta_250_450_' + label_name[0:-3]
		stcs2 = mne.read_source_estimate(out_file2)
		stc_subtract=np.subtract(stcs1,stcs2)
		out_file=data_path + 'GrandAverage_SemLoc_Subtract_Coherence_AlphaBeta_250_450_' + label_name[0:-3]
		stc_grand.save(out_file)
"""
