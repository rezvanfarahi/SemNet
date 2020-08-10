"""
=================================================================
Permutation t-test on source data with spatio-temporal clustering
=================================================================

Tests if the evoked response is significantly different between
conditions across subjects (simulated here using one subject's data).
The multiple comparisons problem is addressed with a cluster-level
permutation test across space and time.

"""

# Authors: Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
#          Eric Larson <larson.eric.d@gmail.com>
# License: BSD (3-clause)

print(__doc__)
# Russell's addition
import sys
sys.path.append('/imaging/local/software/python_packages/nibabel/1.3.0')
sys.path.append('/imaging/local/software/python_packages/pysurfer/v0.4')
# End

# for qsub
# (add ! here if needed) /imaging/local/software/anaconda/latest/x86_64/bin/python
sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
###


import os.path as op
import numpy as np
from numpy.random import randn
from scipy import stats as stats

import mne
from mne import (io, spatial_tris_connectivity, compute_morph_matrix,
                 grade_to_tris)
from mne.epochs import equalize_epoch_counts
from mne.stats import (spatio_temporal_cluster_1samp_test,
                       summarize_clusters_stc)
from mne.minimum_norm import apply_inverse, read_inverse_operator
from mne.datasets import sample
from mne.viz import mne_analyze_colormap

###############################################################################
# Set parameters
data_path = '/imaging/rf02/TypLexMEG/'	# where subdirs for MEG data are
inv_path = '/imaging/rf02/TypLexMEG/'
inv_fname = 'InverseOperator_fft_1_48_clean_ica_EMEG-inv.fif'
event_path = '/group/erp/data/olaf.hauk/MEG/TypLexM/data/MF22_new/'
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
n_subjects=17
stim_delay = 0.034 # delay in s ((NOTE! not in the example, taken from Olaf's script. Rezvan))
tmin, tmax = -0.5, 0.7
event_ids = {'cncrt_wrd': 1, 'abs_wrd': 2}
X=np.zeros((20484,35, 17,2))
for ii, meg in enumerate(ll):
	print ii
	raw_fname = data_path + meg + 'semloc_ssstf_fft_1_48_clean_ica_raw.fif' #clean_ssp_
	event_fname = event_path + meg + 'semloc_raw_ssstnew.txt'
	subjects_dir = data_path 

	tmin = -0.5
	tmax = 0.7  # Use a lower tmax to reduce multiple comparisons

	#   Setup for reading the raw data
	raw = io.Raw(raw_fname)
	events = mne.read_events(event_fname)

	###############################################################################
	# Read epochs for all channels, removing a bad one
	print "epochs"
	picks = mne.pick_types(raw.info, eeg=True, meg=True, eog=True, exclude='bads')
	event_id = 1  # concrete
	if ii==3:
    		reject = dict(eeg=150e-6, grad=200e-12, mag=4e-12) #eog=150e-6,
	else:
	    	reject = dict(eeg=120e-6, grad=200e-12, mag=4e-12) #eog=150e-6,
	epochs1 = mne.Epochs(raw, events, event_id, tmin, tmax, picks=picks,
		             baseline=(None, 0), proj=True, reject=reject, preload=True)

	event_id = 2  # abstract
	epochs2 = mne.Epochs(raw, events, event_id, tmin, tmax, picks=picks,
		             baseline=(None, 0), proj=True, reject=reject, preload=True)

	#    Equalize trial counts to eliminate bias (which would otherwise be
	#    introduced by the abs() performed below)
	equalize_epoch_counts([epochs1, epochs2])

	###############################################################################
	print "Transform to source space"

	fname_inv = inv_path + meg + inv_fname
	snr = 3.0
	lambda2 = 1.0 / snr ** 2
	method = "MNE"  # use dSPM method (could also be MNE or sLORETA)
	inverse_operator = read_inverse_operator(fname_inv)
	sample_vertices = [s['vertno'] for s in inverse_operator['src']]

	#    Let's average and compute inverse, resampling to speed things up
	#for cc in range(34):
	#	b1=cc*20+501; b2=(cc+1)*20+500
	#	X[:,cc,ii]=stc_subtract.copy().crop(b1,b2).mean().data.squeeze()
	evoked1 = epochs1.average()
	evoked1.resample(50)
	condition1 = apply_inverse(evoked1, inverse_operator, lambda2, method)
	evoked2 = epochs2.average()
	evoked2.resample(50)
	condition2 = apply_inverse(evoked2, inverse_operator, lambda2, method)

	#    Let's only deal with t > 0, cropping to reduce multiple comparisons
	condition1.crop(0, None)
	condition2.crop(0, None)
	tmin = condition1.tmin
	tstep = condition1.tstep

	###############################################################################
	print "transform to common cortical space"

	#    Normally you would read in estimates across several subjects and morph
	#    them to the same cortical space (e.g. fsaverage). For example purposes,
	#    we will simulate this by just having each "subject" have the same
	#    response (just noisy in source space) here. Note that for 7 subjects
	#    with a two-sided statistical test, the minimum significance under a
	#    permutation test is only p = 1/(2 ** 6) = 0.015, which is large.
	n_vertices_sample, n_times = condition1.data.shape
	X1 = randn(n_vertices_sample, n_times, 1, 2)
	X1[:, :, :, 0]= condition1.data[:, :, np.newaxis]
	X1[:, :, :, 1]= condition2.data[:, :, np.newaxis] 
	fsave_vertices = [np.arange(10242), np.arange(10242)]
	morph_mat = compute_morph_matrix(subjects[ii], 'avgsubj', sample_vertices,
                                 fsave_vertices, 20, subjects_dir)
	n_vertices_fsave = morph_mat.shape[0]

	#    We have to change the shape for the dot() to work properly
	X1 = X1.reshape(n_vertices_sample, n_times * 1 * 2)
	print('Morphing data.')
	X1 = morph_mat.dot(X1)  # morph_mat is a sparse matrix
	X1 = X1.reshape(n_vertices_fsave, n_times, 1, 2)
	X[:,:,ii,:]=X1[:,:,0,:]

	

	#    It's a good idea to spatially smooth the data, and for visualization
	#    purposes, let's morph these to fsaverage, which is a grade 5 source space
	#    with vertices 0:10242 for each hemisphere. Usually you'd have to morph
	#    each subject's data separately (and you might want to use morph_data
	#    instead), but here since all estimates are on 'sample' we can use one
	#    morph matrix for all the heavy lifting.

	#    Finally, we want to compare the overall activity levels in each condition,
	#    the diff is taken along the last axis (condition). The negative sign makes
	#    it so condition1 > condition2 shows up as "red blobs" (instead of blue).
X2 = np.abs(X)  # only magnitude
X = X[:, :, :, 0] - X[:, :, :, 1]  # make paired contrast
X2 = X2[:, :, :, 0] - X2[:, :, :, 1]  # make paired contrast


	###############################################################################
# Compute statistic

#    To use an algorithm optimized for spatio-temporal clustering, we
#    just pass the spatial connectivity matrix (instead of spatio-temporal)
print('Computing connectivity.')
connectivity = spatial_tris_connectivity(grade_to_tris(5))

#    Note that X needs to be a multi-dimensional array of shape
#    samples (subjects) x time x space, so we permute dimensions
X = np.transpose(X, [2, 1, 0])
X2 = np.transpose(X2, [2, 1, 0])

#    Now let's actually do the clustering. This can take a long time...
#    Here we set the threshold quite high to reduce computation.
p_threshold = 0.05
t_threshold = -stats.distributions.t.ppf(p_threshold / 2., n_subjects - 1)
print('Clustering.')
T_obs, clusters, cluster_p_values, H0 = clu = \
    spatio_temporal_cluster_1samp_test(X, connectivity=connectivity, n_jobs=2,
	                               threshold=t_threshold)
#    Now select the clusters that are sig. at p < 0.05 (note that this value
#    is multiple-comparisons corrected).
good_cluster_inds = np.where(cluster_p_values < 0.1)[0]

###############################################################################
# Visualize the clusters

print('Visualizing clusters.')

#    Now let's build a convenient representation of each cluster, where each
#    cluster becomes a "time point" in the SourceEstimate
stc_all_cluster_vis = summarize_clusters_stc(clu, tstep=tstep,
	                                     vertno=fsave_vertices,
	                                     subject='avgsubj')

out_file1=data_path + 'Permutation_clusterp05_p05_17subj_SemLoc_icaclean_ConcreteAbstract_0_600_noneqlizedEvokednew'
stc_all_cluster_vis.save(out_file1)

Matx=np.zeros((20484,35))
Matx[clusters[good_cluster_inds][1],clusters[good_cluster_inds][0]]=1
tmin1=0
tstep1=20
vertices_to = [np.arange(10242), np.arange(10242)]
matx_stc = mne.SourceEstimate(Matx, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='avgsubj')
out_file2=data_path + 'Permutation_stepwise_clusterp05_p05_17subj_SemLoc_icaclean_ConcreteAbstract_0_600_noneqlizedEvokednew'
plv_stc.save(out_file2)
