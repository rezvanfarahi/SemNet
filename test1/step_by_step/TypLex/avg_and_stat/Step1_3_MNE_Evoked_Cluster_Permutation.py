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
#sys.path.insert(1,'/imaging/rf02/TypLexMEG/mne_python0.9')
#sys.path.insert(1,'/imaging/local/software/mne_python/latest')
sys.path.append('/imaging/local/software/mne_python/latest_v0.9full')

sys.path.insert(1,'/imaging/local/software/python_packages/scikit-learn/v0.13.1/lib.linux-x86_64-2.7/sklearn/externals')
import joblib
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
from mne.stats import permutation_t_test

###############################################################################
# Set parameters
out_path = '/imaging/rf02/TypLexMEG/icaanalysis_results/typlex/stc/permutation/evoked/' # root
data_path = '/imaging/rf02/TypLexMEG/'	# where subdirs for MEG data are
inv_path = '/imaging/rf02/TypLexMEG/'
inv_fname = 'typlex_InvOp_ico5diagreg_fft_1_48_clean_ica_EMEG-inv.fif'#inv_fname = 'InverseOperator_fft_1_48_clean_ica_EMEG-inv.fif'
event_path = '/group/erp/data/olaf.hauk/MEG/TypLexM/data/MF22_new/'
label_path = '/imaging/rf02/TypLexMEG/fsaverage/label'

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
X=np.zeros((20484,5, 17))
for p_threshold in ll:
    for ii, meg in enumerate(list_all):
    	print ii
    	fname_cncrt = data_path + meg + 'firstMorphed_ico_typlex_ica_typword_Source_Evoked_m500_700'
    	stc_cncrt = mne.read_source_estimate(fname_cncrt)
    	#abslt_cncrt= np.absolute(stc_cncrt.copy().crop(.3,.4).data)
    	fname_abs = data_path + meg + 'firstMorphed_ico_typlex_ica_typnonword_Source_Evoked_m500_700'
    	stc_abs = mne.read_source_estimate(fname_abs)
    	#    Let's only deal with t > 0, cropping to reduce multiple comparisons
#    	stc_cncrt.resample(100)
#    	print "C fin"
#    	stc_abs.resample(100)
#    	print "A fin"
    	stc_subtract=np.subtract(stc_cncrt,stc_abs)

    	stc_subtract_data=np.subtract(np.abs(stc_cncrt.data),np.abs(stc_abs.data))
    	b2=0.05
    	for cc in range(5):
    		b1=b2+0.0001; b2=b1+0.1-0.0001; print b1; print b2
    	 	X[:,cc,ii]=np.subtract(np.abs(stc_cncrt.copy().crop(b1,b2).mean().data.squeeze()),np.abs(stc_abs.copy().crop(b1,b2).mean().data.squeeze()))#stc_subtract.copy().crop(b1,b2).mean().data.squeeze()
    	#X[:,:,ii]=stc_subtract_data

#    	tmin = stc_subtract.tmin
#     
#    	tstep = stc_subtract.tstep
    	#X[:,1,:]=X[:,0,:].copy()
    	###############################################################################
    	print "transform to common cortical space"
    
    
    print('Computing connectivity.')
    fname_label = label_path + '/' + 'toremove_wbspokes-lh.label'; labelL = mne.read_label(fname_label)
    fname_label = label_path + '/' + 'toremove_wbspokes-rh.label'; labelR = mne.read_label(fname_label)
    labelss=labelL+labelR
    bb=stc_cncrt.in_label(labelss)
    fsave_vertices = [np.arange(10242), np.arange(10242)]
    nnl=np.in1d(fsave_vertices[0],bb.lh_vertno)
    nnr=np.in1d(fsave_vertices[1],bb.rh_vertno)
    spatial_exclude=np.hstack((fsave_vertices[0][nnl], fsave_vertices[0][nnr]+10242))
    
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
#    def stat_fun(*args):
#        return permutation_t_test
    print('Clustering.')
    tail=0
    max_step=5;
    T_obs, clusters, cluster_p_values, H0 = clu = spatio_temporal_cluster_1samp_test(X, connectivity=connectivity,  threshold=t_threshold, tail=tail, t_power=1, step_down_p=0.05,spatial_exclude=spatial_exclude)#,stat_fun=permutation_t_test)#, max_step=0)#, max_step=0)#spatial_exclude=exclude_ver, 

    #    Now select the clusters that are sig. at p < 0.05 (note that this value
    #    is multiple-comparisons corrected).
    p_thresh=cluster_p_values.min()+0.0001#0.05
    good_cluster_inds = np.where(cluster_p_values <p_thresh)[0]
    print cluster_p_values[good_cluster_inds]; print good_cluster_inds
    ###############################################################################
    # Visualize the clusters
    
    print('Visualizing clusters.')
    
    #    Now let's build a convenient representation of each cluster, where each
    #    cluster becomes a "time point" in the SourceEstimate
    stc_all_cluster_vis = summarize_clusters_stc(clu, tstep=0.1, vertices=fsave_vertices, subject='fsaverage', p_thresh=p_thresh)
    
    out_file1=out_path + 'ClusPer_abs_morphedico_cp'+str(p_threshold)[2:]+'_p'+str(p_thresh)[2:]+'_17subj_typlex_ica_typWordPW_50_550_100ms_tail'+str(tail)+'_maxstep'+str(max_step)
    stc_all_cluster_vis.save(out_file1)
    
    Matx=np.zeros((20484,5))
    T=np.divide(T_obs,np.absolute(T_obs))
    for cc in range(good_cluster_inds.shape[0]):
    	Matx[clusters[good_cluster_inds[cc]][1],clusters[good_cluster_inds[cc]][0]]=T[clusters[good_cluster_inds[cc]][0],clusters[good_cluster_inds[cc]][1]]
    	
    
    #aa=Matx[clusters[good_cluster_inds[cc]][1],clusters[good_cluster_inds[cc]][0]]
    #bb=T_obs[clusters[good_cluster_inds[cc]][0],clusters[good_cluster_inds[cc]][1]]<0
    #aa[bb]=-1
    #Matx[clusters[good_cluster_inds[cc]][1],clusters[good_cluster_inds[cc]][0]]=aa
    
    tmin1=50
    tstep1=100
    vertices_to = [np.arange(10242), np.arange(10242)]
    matx_stc = mne.SourceEstimate(Matx, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
    out_file2=out_path + 'ClusPer_abs_sw_morphedico_cp'+str(p_threshold)[2:]+'_p'+str(p_thresh)[2:]+'_17subj_typlex_ica_typWordPW_50_550_100ms_tail'+str(tail)+'_maxstep'+str(max_step)
    matx_stc.save(out_file2)

