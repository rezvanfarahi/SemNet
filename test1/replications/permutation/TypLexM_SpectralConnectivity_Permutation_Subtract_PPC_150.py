
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
#sys.path.append('/imaging/local/software/python_packages/nibabel/1.3.0')
sys.path.insert(1,'/imaging/rf02/scikit-learn-0.16.1')
#sys.path.append('/imaging/local/software/python_packages/scikit-learn/v0.14.1')
sys.path.append('/imaging/local/software/mne_python/latest_v0.9')
# for qsub
# (add ! here if needed) 
#sys.path.insert(1,'/imaging/local/software/anaconda/1.9.1/x86_64/envs/py3k')#bin/python')
sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
###
import numpy as np
import mne
from mne import io

import os
import scipy.io

from mne import (io, spatial_tris_connectivity, compute_morph_matrix,
                 grade_to_tris)
from mne.stats import (spatio_temporal_cluster_1samp_test,
                       summarize_clusters_stc)


###############################################################################
data_path = '/imaging/rf02/TypLexMEG/' # root directory for your MEG data
os.chdir(data_path)
subjects_dir = '/home/rf02/rezvan/TypLexMEG/'    # where your MRI subdirectories are
# where event-files are
out_path = '/imaging/rf02/TypLexMEG/icaanalysis_results/stc/permutation/connectivity/WB_spokes/'    # where event files are

label_path = '/home/rf02/rezvan/TypLexMEG/createdlabels_SL/'

inv_path = '/imaging/rf02/TypLexMEG/'
inv_fname = 'InverseOperator_EMEG-inv.fif'

# get indices for subjects to be processed from command line input
# 
print sys.argv
p_inds = []
for ss in sys.argv[1:]:
   p_inds.append( int( ss ) )

subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16] # removed
#p_inds=[0,1,2,3,4,5,6,7,8,9,10]
p_list=[0.001,0.05,0.045,0.04,0.03,0.025,0.02,0.01,0.008,0.005,0.002,0.001,0.0005,0.0001]
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
for ss in p_inds:
 ll.append(p_list[ss])

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
labellist = [ 'lh.precentral','rh.precentral','lh.postcentral','rh.postcentral', 'lh.lateraloccipital','rh.lateraloccipital']#labellist = [ 'posterior_fusiform-lh','rh.lateralorbitofron-rh', 'rh.medialorbitofron-rh']#labellist = ['leftATL-rh','rightATL-lh']#, 'rightATL-rh']#'atlleft-lh']#, 'atlright-rh']#,'medialtempright-rh','medialtempleft-lh']
trange='_150_350_'
conn_method='PPC'
for p_threshold in ll:
    for label_name in labellist:
    	ii=-1
    	#X=np.zeros((20484,2,17))
    	X_theta=np.zeros((20484,17))
    	X_alpha=np.zeros((20484,17))
    	X_beta=np.zeros((20484,17))
    	X_gamma=np.zeros((20484,17))
    	for meg in list_all:
    		ii=ii+1;
    		#print ii
    	
    ### Concrete
    		fname_cncrt = data_path + meg + 'firstMorphed_SemLoc_icaclean_equalized_Concrete_' + conn_method + '_ThetaAlphaBetaGamma' + trange + label_name[0:-3] 
    		stc_cncrt = mne.read_source_estimate(fname_cncrt)
    		fname_abs = data_path + meg + 'firstMorphed_SemLoc_icaclean_equalized_Abstract_' + conn_method + '_ThetaAlphaBetaGamma' + trange + label_name[0:-3] 
    		stc_abs = mne.read_source_estimate(fname_abs)	
    		stc_data_subtract=np.subtract(np.abs(stc_cncrt.data),np.abs(stc_abs.data))
    		#X[:,:,ii]=stc_subtract.data
    		X_theta[:,ii]=stc_data_subtract[:,0]
    		X_alpha[:,ii]=stc_data_subtract[:,1]
    		X_beta[:,ii]=stc_data_subtract[:,2]
    		X_gamma[:,ii]=stc_data_subtract[:,3]
    		#X[:,:,ii]=stc_subtract.data	
    		
    	#X = np.transpose(X, [2, 1, 0])	
    	connectivity = spatial_tris_connectivity(grade_to_tris(5))
    	#p_threshold = 0.045
    	print p_threshold
    	n_subjects=17
    	t_threshold = -scipy.stats.distributions.t.ppf(p_threshold/2., n_subjects - 1)
    	fsave_vertices = [np.arange(10242), np.arange(10242)]
    	n_permutations=1024
    	# I'm replicating the columns because it makes the script faster. Have checked and it doesn't introduce 
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
	print('Computing connectivity.')
	labelss = mne.read_labels_from_annot(subject='fsaverage', parc='aparc', subjects_dir=data_path)
	labelss.remove(labelss[-2])
	labelss.remove(labelss[-1])
	bb=stc_cncrt.in_label(np.sum(labelss))
	fsave_vertices = [np.arange(10242), np.arange(10242)]
	nnl=~np.in1d(fsave_vertices[0],bb.lh_vertno)
	nnr=~np.in1d(fsave_vertices[1],bb.rh_vertno)
	spatial_exclude=np.hstack((fsave_vertices[0][nnl], fsave_vertices[0][nnr]+10242))
        
    
    	
    	print('Clustering theta')
    	T_obst, clusterst, cluster_p_valuest, H0t = clut = spatio_temporal_cluster_1samp_test(X_theta, connectivity=connectivity,  n_permutations=n_permutations, n_jobs=4, threshold=t_threshold, tail=0,t_power=1, spatial_exclude=spatial_exclude,step_down_p=0.05)
    	#    Now select the clusters that are sig. at p < 0.05 (note that this value
    	#    is multiple-comparisons corrected).
    	good_cluster_indst = np.where(cluster_p_valuest < 0.05)[0]
    	print cluster_p_valuest.min()
    	
    	
    	print('Clustering alpha')
    	T_obsa, clustersa, cluster_p_valuesa, H0a = clua = spatio_temporal_cluster_1samp_test(X_alpha, connectivity=connectivity,  n_permutations=n_permutations, n_jobs=4, threshold=t_threshold, tail=0,t_power=1, spatial_exclude=spatial_exclude,step_down_p=0.05)
    	#    Now select the clusters that are sig. at p < 0.05 (note that this value
    	#    is multiple-comparisons corrected).
    	good_cluster_indsa = np.where(cluster_p_valuesa < 0.05)[0]
    	print cluster_p_valuesa.min()
    	
    	print('Clustering beta')
    	T_obsb, clustersb, cluster_p_valuesb, H0b = club = spatio_temporal_cluster_1samp_test(X_beta, connectivity=connectivity,  n_permutations=n_permutations, n_jobs=4, threshold=t_threshold, tail=0,t_power=1,spatial_exclude=spatial_exclude, step_down_p=0.05)
    	#    is multiple-comparisons corrected).
    	good_cluster_indsb = np.where(cluster_p_valuesb < 0.05)[0]
    	print cluster_p_valuesb.min()
    	
    	print('Clustering gamma')
    	T_obsg, clustersg, cluster_p_valuesg, H0g = clug = spatio_temporal_cluster_1samp_test(X_gamma, connectivity=connectivity,  n_permutations=n_permutations, n_jobs=4, threshold=t_threshold, tail=0,t_power=1,spatial_exclude=spatial_exclude, step_down_p=0.05)
    	#    Now select the clusters that are sig. at p < 0.05 (note that this value
    	#    is multiple-comparisons corrected).
    	good_cluster_indsg = np.where(cluster_p_valuesg < 0.05)[0]
    	print cluster_p_valuesg.min()
    	
    
    	tstep=1e-3
    	fsave_vertices = [np.arange(10242), np.arange(10242)]
    	p_thresh=0.05#cluster_p_valuest.min()+0.000001
    	
    	if cluster_p_valuest.min()<=0.05:
    		print p_thresh

    		stc_all_cluster_vist = summarize_clusters_stc(clut, tstep=tstep, p_thresh=p_thresh, vertno=fsave_vertices, subject='fsaverage')
    		out_file1=out_path + 'Permutation_firstMorphed_SemLoc_icaclean_Subtract_' + conn_method + '_Theta' + trange + '_' + str(p_threshold)[2:5] + '_' + str(p_thresh)[2:5]+ '_' + label_name[0:-3]
    		stc_all_cluster_vist.save(out_file1)
    
    	p_thresh=0.05#cluster_p_valuesa.min()+0.000001
    	
    	if cluster_p_valuesa.min()<=0.05:
    		print p_thresh
    		stc_all_cluster_visa = summarize_clusters_stc(clua, tstep=tstep, p_thresh=p_thresh, vertno=fsave_vertices, subject='fsaverage')
    		out_file1=out_path + 'Permutation_firstMorphed_SemLoc_icaclean_Subtract_' + conn_method + '_Alpha' + trange + '_' + str(p_threshold)[2:5]+ '_' + str(p_thresh)[2:5] + '_' + label_name[0:-3]
    		stc_all_cluster_visa.save(out_file1)
    
    	p_thresh=0.05#cluster_p_valuesb.min()+0.000001
    	if cluster_p_valuesb.min()<=0.05:
    		print p_thresh
    		stc_all_cluster_visb = summarize_clusters_stc(club, tstep=tstep,p_thresh=p_thresh, vertno=fsave_vertices, subject='fsaverage')
    		out_file1=out_path + 'Permutation_firstMorphed_SemLoc_icaclean_Subtract_' + conn_method + '_Beta' + trange + '_' + str(p_threshold)[2:5]+ '_' + str(p_thresh)[2:5] + '_' + label_name[0:-3]
    		stc_all_cluster_visb.save(out_file1)
    
    	p_thresh=0.05#cluster_p_valuesg.min()+0.000001
    	if cluster_p_valuesg.min()<=0.05:
    		print p_thresh
    		stc_all_cluster_visg = summarize_clusters_stc(clug, tstep=tstep, p_thresh=p_thresh, vertno=fsave_vertices, subject='fsaverage')    
    		out_file1=out_path + 'Permutation_firstMorphed_SemLoc_icaclean_Subtract_' + conn_method + '_Gamma' + trange + '_' + str(p_threshold)[2:5]+ '_' + str(p_thresh)[2:5] + '_' + label_name[0:-3]
    		stc_all_cluster_visg.save(out_file1)
    
