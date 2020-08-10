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
sys.path.insert(1,'/imaging/local/software/mne_python/latest_v0.9')

import os.path as op
import numpy as np
from numpy.random import randn
from scipy import stats as stats

import mne
from mne import (io, spatial_tris_connectivity, compute_morph_matrix,spatio_temporal_tris_connectivity,
                 grade_to_tris)
from mne.epochs import equalize_epoch_counts
from mne.stats import (spatio_temporal_cluster_1samp_test,
                       summarize_clusters_stc)
from mne.minimum_norm import apply_inverse, read_inverse_operator
from mne.datasets import sample
from mne.viz import mne_analyze_colormap

###############################################################################
# Set parameters
out_path = '/imaging/rf02/TypLexMEG/icaanalysis_results/stc/permutation/evoked/' # root
data_path = '/imaging/rf02/TypLexMEG/'	# where subdirs for MEG data are
inv_path = '/imaging/rf02/TypLexMEG/'
inv_fname = 'InverseOperator_fft_1_48_clean_ica_EMEG-inv.fif'
event_path = '/group/erp/data/olaf.hauk/MEG/TypLexM/data/MF22_new/'
print sys.argv
p_inds = []
for ss in sys.argv[1:]:
   p_inds.append( int( ss ) )

subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16] # removed
#p_inds=[0,1,2,3,4,5,6,7,8,9,10]
p_list=[0.05,0.045,0.04,0.03,0.025,0.02,0.01,0.008,0.005,0.002,0.001, 0.0008,0.0005,0.0002,0.0001,0.00005,0.00001]#0.001,0.0009,0.0008,0.0007,0.0006,0.0005]#
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

n_subjects=17
stim_delay = 0.034 # delay in s ((NOTE! not in the example, taken from Olaf's script. Rezvan))
tmin, tmax = -0.5, 0.7
event_ids = {'cncrt_wrd': 1, 'abs_wrd': 2}
X=np.zeros((20484,2, 17))
for p_threshold in ll:
    for ii, meg in enumerate(list_all):
    	print ii
    	fname_cncrt = data_path + meg + 'firstMorphed_SemLoc_icaclean_Concrete_Source_Evoked_m500_700'
    	stc_cncrt = mne.read_source_estimate(fname_cncrt)
    	#abslt_cncrt= np.absolute(stc_cncrt.copy().crop(.3,.4).data)
    	fname_abs = data_path + meg + 'firstMorphed_SemLoc_icaclean_Abstract_Source_Evoked_m500_700'
    	stc_abs = mne.read_source_estimate(fname_abs)
    	#    Let's only deal with t > 0, cropping to reduce multiple comparisons
    	stc_subtract=np.subtract(stc_cncrt,stc_abs)
    	b2=0.2
    	for cc in range(2):
    		b1=b2+0.0001; b2=b1+0.1-0.0001; print b1; print b2
    	 	X[:,cc,ii]=np.subtract(np.abs(stc_cncrt.copy().crop(b1,b2).mean().data.squeeze()),np.abs(stc_abs.copy().crop(b1,b2).mean().data.squeeze()))#stc_subtract.copy().crop(b1,b2).mean().data.squeeze()
    	tmin = stc_subtract.tmin
     
    	tstep = stc_subtract.tstep
    	X[:,1,:]=X[:,0,:].copy()
    	###############################################################################
    	print "transform to common cortical space"
    
    	#    Normally you would read in estimates across several subjects and morph
    	#    them to the same cortical space (e.g. fsaverage). For example purposes,
    	#    we will simulate this by just having each "subject" have the same
    	#    response (just noisy in source space) here. Note that for 7 subjects
    	#    with a two-sided statistical test, the minimum significance under a
    	#    permutation test is only p = 1/(2 ** 6) = 0.015, which is large.
    	
    	#    It's a good idea to spatially smooth the data, and for visualization
    	#    purposes, let's morph these to fsaverage, which is a grade 5 source space
    	#    with vertices 0:10242 for each hemisphere. Usually you'd have to morph
    	#    each subject's data separately (and you might want to use morph_data
    	#    instead), but here since all estimates are on 'sample' we can use one
    	#    morph matrix for all the heavy lifting.
    
    	#    Finally, we want to compare the overall activity levels in each condition,
    	#    the diff is taken along the last axis (condition). The negative sign makes
    	#    it so condition1 > condition2 shows up as "red blobs" (instead of blue).
    print('Computing connectivity.')
    labelss = mne.read_labels_from_annot(subject='fsaverage', parc='rezvan_annot', subjects_dir=data_path)
    labelss.remove(labelss[-2])
    labelss.remove(labelss[-1])
    bb=stc_cncrt.in_label(np.sum(labelss))
    fsave_vertices = [np.arange(10242), np.arange(10242)]
    nnl=~np.in1d(fsave_vertices[0],bb.lh_vertno)
    nnr=~np.in1d(fsave_vertices[1],bb.rh_vertno)
    exclude_ver=np.hstack((fsave_vertices[0][nnl], fsave_vertices[0][nnr]+10242))
    
    #to_keep=np.hstack((bb.lh_vertno, bb.rh_vertno+10242))
    #gts=grade_to_tris(5)
    #to_con=np.zeros((1,3))
    #for cntr in range(gts.shape[0]):
    #    if all(np.in1d(gts[cntr,:],to_keep)):
    #        to_con=np.vstack((to_con,gts[cntr,:]))
    #to_con=np.delete(to_con, (0), axis=0)
    connectivity = spatial_tris_connectivity(grade_to_tris(5))   
    #connectivity = spatio_temporal_tris_connectivity(grade_to_tris(5), n_times=5)   
    n_permutations=1024
    
    #    Note that X needs to be a multi-dimensional array of shape
    #    samples (subjects) x time x space, so we permute dimensions
    X = np.transpose(X, [2, 1, 0])
    
    #    Now let's actually do the clustering. This can take a long time...
    #    Here we set the threshold quite high to reduce computation.
    cstep=0.1#(np.abs(X).max()-np.abs(X).min())/100
    cstart=0#np.abs(X).min()
    #p_threshold = 0.043
    t_threshold = -stats.distributions.t.ppf(p_threshold/2., n_subjects - 1)#dict(start=0, step=.1)#
    #t_threshold=2
    print('Clustering.')
    tail=0
    T_obs, clusters, cluster_p_values, H0 = clu = spatio_temporal_cluster_1samp_test(X, connectivity=connectivity,  threshold=t_threshold, tail=tail, t_power=1, step_down_p=0.05)#, max_step=0)#, max_step=0)#spatial_exclude=exclude_ver, 
    #    Now select the clusters that are sig. at p < 0.05 (note that this value
    #    is multiple-comparisons corrected).
    p_thresh=0.05
    good_cluster_inds = np.where(cluster_p_values <p_thresh)[0]
    print cluster_p_values.min(); print good_cluster_inds
    ###############################################################################
    # Visualize the clusters
    
    print('Visualizing clusters.')
    
    #    Now let's build a convenient representation of each cluster, where each
    #    cluster becomes a "time point" in the SourceEstimate
    stc_all_cluster_vis = summarize_clusters_stc(clu, tstep=0.1, vertno=fsave_vertices, subject='fsaverage', p_thresh=p_thresh)
    
    out_file1=out_path + 'ClusterPermutation_abs_firstmorphed_clusterp'+str(p_threshold)[2:]+'_p'+str(p_thresh)[2:]+'_17subj_SemLoc_icaclean_ConcreteAbstract_200_300_tail'+str(tail)
    stc_all_cluster_vis.save(out_file1)
    
#    Matx=np.zeros((20484,2))
#    T=np.divide(T_obs,np.absolute(T_obs))
#    for cc in range(good_cluster_inds.shape[0]):
#    	Matx[clusters[good_cluster_inds[cc]][1],clusters[good_cluster_inds[cc]][0]]=T[clusters[good_cluster_inds[cc]][0],clusters[good_cluster_inds[cc]][1]]
#    	
#    
#    #aa=Matx[clusters[good_cluster_inds[cc]][1],clusters[good_cluster_inds[cc]][0]]
#    #bb=T_obs[clusters[good_cluster_inds[cc]][0],clusters[good_cluster_inds[cc]][1]]<0
#    #aa[bb]=-1
#    #Matx[clusters[good_cluster_inds[cc]][1],clusters[good_cluster_inds[cc]][0]]=aa
#    
#    tmin1=0
#    tstep1=100
#    vertices_to = [np.arange(10242), np.arange(10242)]
#    matx_stc = mne.SourceEstimate(Matx, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
#    out_file2=out_path + 'ClusterPermutation_abs_stepwise_firstmorphed_clusterp'+str(p_threshold)[2:]+'_p'+str(p_thresh)[2:]+'_17subj_SemLoc_icaclean_ConcreteAbstract_200_300_100ms_tail'+str(tail)
#    matx_stc.save(out_file2)

