
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
from mne.layouts import read_layout
from mne.time_frequency import compute_raw_psd
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

subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16] # removed
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
reject = dict(eeg=120e-6, grad=200e-12, mag=4e-12) #eog=150e-6, 
event_ids = {'cncrt_wrd': 1, 'abs_wrd': 2}
# frequency bands with variable number of cycles for wavelets
stim_delay = 0.034 # delay in s
frequencies = np.arange(8, 45, 1)
n_cycles = frequencies / float(7)
n_cycles[frequencies<=20] = 2
n_cycles[frequencies<=20] += np.arange(0,1,0.077)
    # n_cycles[np.logical_and((frequencies>10), (frequencies<=20))] = 3
labellist = ['atlleft-lh']#, 'atlright-rh']#,'medialtempright-rh','medialtempleft-lh']


for label_name in labellist:
	ii=-1
	#X=np.zeros((20484,2,17))
	X_theta=np.zeros((20484,17))
	X_alpha=np.zeros((20484,17))
	X_beta=np.zeros((20484,17))
	X_gamma=np.zeros((20484,17))
	for meg in ll:
		ii=ii+1;
		#print ii
	
### Concrete
		fname_cncrt = data_path + meg + 'Morphed_SemLoc_icaclean_Concrete_PhaseLocking_ThetaAlphaBetaGamma_50_250_' + label_name[0:-3]
		stc_cncrt = mne.read_source_estimate(fname_cncrt)
		fname_abs = data_path + meg + 'Morphed_SemLoc_icaclean_Abstract_PhaseLocking_ThetaAlphaBetaGamma_50_250_' + label_name[0:-3]
		stc_abs = mne.read_source_estimate(fname_abs)	
		stc_subtract=np.subtract(stc_cncrt,stc_abs)
		#X[:,:,ii]=stc_subtract.data
		X_theta[:,ii]=stc_subtract.data[:,0]
		X_alpha[:,ii]=stc_subtract.data[:,1]
		X_beta[:,ii]=stc_subtract.data[:,2]
		X_gamma[:,ii]=stc_subtract.data[:,3]	
		
	#X = np.transpose(X, [2, 1, 0])	
	connectivity = spatial_tris_connectivity(grade_to_tris(5))
	p_threshold = 0.0008
	print p_threshold
	n_subjects=17
	t_threshold = -scipy.stats.distributions.t.ppf(p_threshold / 2., n_subjects - 1)
	fsave_vertices = [np.arange(10242), np.arange(10242)]
	n_permutations=1024

	xt=np.zeros((17,2,20484))
	xt[:,0,:]=X_theta.transpose()
	xt[:,1,:]=X_theta.transpose()
	X_theta=xt

	xa=np.zeros((17,2,20484))
	xa[:,0,:]=X_alpha.transpose()
	xa[:,1,:]=X_alpha.transpose()
	X_alpha=xa

	xb=np.zeros((17,2,20484))
	xb[:,0,:]=X_beta.transpose()
	xb[:,1,:]=X_beta.transpose()
	X_beta=xb

	xg=np.zeros((17,2,20484))
	xg[:,0,:]=X_gamma.transpose()
	xg[:,1,:]=X_gamma.transpose()
	X_gamma=xg


	
	print('Clustering theta')
	T_obst, clusterst, cluster_p_valuest, H0t = clut = spatio_temporal_cluster_1samp_test(X_theta, connectivity=connectivity, n_jobs=1, threshold=t_threshold)
	#    Now select the clusters that are sig. at p < 0.05 (note that this value
	#    is multiple-comparisons corrected).
	good_cluster_indst = np.where(cluster_p_valuest < 0.05)[0]
	print cluster_p_valuest.min()
	

	print('Clustering alpha')
	T_obsa, clustersa, cluster_p_valuesa, H0a = clua = spatio_temporal_cluster_1samp_test(X_alpha, connectivity=connectivity, n_jobs=1, threshold=t_threshold)
	#    Now select the clusters that are sig. at p < 0.05 (note that this value
	#    is multiple-comparisons corrected).
	good_cluster_indsa = np.where(cluster_p_valuesa < 0.05)[0]
	print cluster_p_valuesa.min()

	print('Clustering beta')
	T_obsb, clustersb, cluster_p_valuesb, H0b = club = spatio_temporal_cluster_1samp_test(X_beta, connectivity=connectivity, n_jobs=1, threshold=t_threshold)
	#    Now select the clusters that are sig. at p < 0.05 (note that this value
	#    is multiple-comparisons corrected).
	good_cluster_indsb = np.where(cluster_p_valuesb < 0.05)[0]
	print cluster_p_valuesb.min()

	print('Clustering gamma')
	T_obsg, clustersg, cluster_p_valuesg, H0g = clug = spatio_temporal_cluster_1samp_test(X_gamma, connectivity=connectivity, n_jobs=1, threshold=t_threshold)
	#    Now select the clusters that are sig. at p < 0.05 (note that this value
	#    is multiple-comparisons corrected).
	good_cluster_indsg = np.where(cluster_p_valuesg < 0.05)[0]
	print cluster_p_valuesg.min()
	

	tstep=1e-3
	fsave_vertices = [np.arange(10242), np.arange(10242)]
	p_thresh=cluster_p_valuest.min()+0.001
	print p_thresh
	stc_all_cluster_vist = summarize_clusters_stc(clut, tstep=tstep, p_thresh=p_thresh, vertno=fsave_vertices, subject='avgsubj')
	out_file1=data_path + 'Permutation_SemLoc_icaclean_Concrete_PhaseLocking_Theta_50_250_' + label_name[0:-3]
	stc_all_cluster_vist.save(out_file1)

	p_thresh=cluster_p_valuesa.min()+0.001
	print p_thresh
	stc_all_cluster_visa = summarize_clusters_stc(clua, tstep=tstep, p_thresh=p_thresh, vertno=fsave_vertices, subject='avgsubj')
	out_file1=data_path + 'Permutation_SemLoc_icaclean_Concrete_PhaseLocking_Alpha_50_250_' + label_name[0:-3]
	stc_all_cluster_visa.save(out_file1)

	p_thresh=cluster_p_valuesb.min()+0.001
	print p_thresh
	stc_all_cluster_visb = summarize_clusters_stc(club, tstep=tstep, p_thresh=p_thresh, vertno=fsave_vertices, subject='avgsubj')
	out_file1=data_path + 'Permutation_SemLoc_icaclean_Concrete_PhaseLocking_Beta_50_250_' + label_name[0:-3]
	stc_all_cluster_visb.save(out_file1)

	p_thresh=cluster_p_valuesg.min()+0.001
	print p_thresh
	stc_all_cluster_visg = summarize_clusters_stc(clug, tstep=tstep, p_thresh=p_thresh, vertno=fsave_vertices, subject='avgsubj')

	out_file1=data_path + 'Permutation_SemLoc_icaclean_Concrete_PhaseLocking_Gamma_50_250_' + label_name[0:-3]
	stc_all_cluster_visg.save(out_file1)
