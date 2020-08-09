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

print __doc__

import sys
sys.path.insert(1,'/imaging/local/software/mne_python/latest_v0.11')

#sys.path.insert(1,'/imaging/rf02/Semnet/mne_python0.9')
#sys.path.insert(1,'/imaging/local/software/mne_python/latest')
sys.path.insert(1,'/imaging/local/software/python_packages/scikit-learn/v0.13.1/lib.linux-x86_64-2.7/sklearn/externals')
import joblib


import os.path as op
import numpy as np
from numpy.random import randn
from scipy import stats as stats

import mne
from mne import (io, spatial_tris_connectivity, compute_morph_matrix,spatio_temporal_tris_connectivity,read_source_estimate,grade_to_tris)
from mne.epochs import equalize_epoch_counts
from mne.stats import (spatio_temporal_cluster_1samp_test,
                       summarize_clusters_stc)
from mne.minimum_norm import apply_inverse, read_inverse_operator
from mne.stats import f_threshold_mway_rm, f_mway_rm, fdr_correction, spatio_temporal_cluster_test
from mne.datasets import sample
from mne.viz import mne_analyze_colormap
import os

###############################################################################
# Set parameters

exclude_wbmedial=False
exclude_ROIs=True
if exclude_wbmedial:
    out_path = '/imaging/rf02/Semnet/semnet4semloc/stc/permutation/evoked/LD/' # root LD or SD
uvttest_path = '/imaging/rf02/Semnet/semnet4semloc/stc/uvttest/evoked/LD/1_48/' # root

if exclude_ROIs:
    out_path = '/imaging/rf02/Semnet/semnet4semloc/stc/permutation/masked_ROIs/evoked/LD/' # root LD or SD

if not os.path.exists(out_path):
    os.makedirs(out_path)
if not os.path.exists(uvttest_path):
    os.makedirs(uvttest_path)
data_path = '/imaging/rf02/Semnet/'	# where subdirs for MEG data are
inv_path = '/imaging/rf02/Semnet/'
print sys.argv
p_inds = []
for ss in sys.argv[1:]:
   p_inds.append( int( ss ) )
label_path = '/imaging/rf02/TypLexMEG/fsaverage/label'

#subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19] # removed
#p_inds=[0,1,2,3,4,5,6,7,8,9,10]
p_list=[0.1,0.05,0.045,0.04,0.025,0.01,0.005]#0.001,0.0009,0.0008,0.0007,0.0006,0.0005]#
#print "subject_inds:"
#print subject_inds
print "No rejection"

list_all =  ['/meg16_0032/160218/', #1 '/meg16_0030/160216/', #0 not needed for LD
            '/meg16_0034/160219/', #3
            '/meg16_0035/160222/', #4
            '/meg16_0042/160229/', #7
            '/meg16_0045/160303/', #8
            '/meg16_0052/160310/', #10
            '/meg16_0056/160314/',#11
            '/meg16_0069/160405/',#12 
            '/meg16_0070/160407/', #13
            '/meg16_0072/160408/', #15
            '/meg16_0073/160411/', #16
            '/meg16_0075/160411/', #17
            '/meg16_0078/160414/', #18
            '/meg16_0082/160418/', #19
            '/meg16_0086/160422/', #20
            '/meg16_0097/160512/', #21 
            '/meg16_0122/160707/', #22 
            '/meg16_0125/160712/', #24 
            ]



# subjects names used for MRI data
subjects=  ['MRI_meg16_0032' ,#1 'MRI_meg16_0030' ,#0 not needed for LD
            'MRI_meg16_0034' ,#2
            'MRI_meg16_0035' ,#3
            'MRI_meg16_0042' ,#4
            'MRI_meg16_0045' ,#5
            'MRI_meg16_0052' ,#6
            'MRI_meg16_0056' ,#7
    	       'MRI_meg16_0069' ,#8
 	       'MRI_meg16_0070' ,#9
            'MRI_meg16_0072' ,#10
            'MRI_meg16_0073' ,#11
            'MRI_meg16_0075' ,#12
            'MRI_meg16_0078' ,#13
            'MRI_meg16_0082' ,#14
            'MRI_meg16_0086' ,#15
            'MRI_meg16_0097' ,#16
            'MRI_meg16_0122' ,#17
            'MRI_meg16_0125' ,#18
            ]

ll = []
for ss in p_inds:
 ll.append(p_list[ss])

print "ll:"
print ll

n_subjects=18
stim_delay = 0.034 # delay in s ((NOTE! not in the example, taken from Olaf's script. Rezvan))
tmin, tmax = -0.3, 0.6
# so we have X: observations * time * vertices in each condition
#X_list=range(2)
event_names = ['Emotional','Concrete']# 'Visual']#'Hear']#, 'Neutral', 'Emotional','Pwordc']'Neutral', 
#thiscond=3
thiscond=1
n_levels=len(event_names)
n_times=4
tmin1=50
tstep1=100
#nwins=5
X=np.zeros((n_subjects,n_times,20484,n_levels))
for p_threshold in ll:
    
    for ii, meg in enumerate(list_all):
        print ii
        for event_no in range(n_levels):     
            fname_in = data_path + meg + 'firstMorphed_ico_oldreg_LD_SL_1_48ica_'+event_names[event_no]+'_Source_Evoked_m300_600'#'firstMorphed_ico_SemDec_pnt1_30ica_'+event_names[event_no]+'_Source_Evoked_mnrsmpl_50_550_100ms_0ov'#'_Source_Evoked_m300_600'#'_Source_Evoked_mnrsmpl_50_550_100ms_0ov'
            stc_cond = read_source_estimate(fname_in)
#            stc_cond.resample(60)
#            stc_cond.crop(0.050,0.450)
#            X[ii,:,:,event_no]=np.transpose(stc_cond.data,[1,0]) 
            wcnt=-1
            for wcnt1,wcnt2 in zip([350,450,550,650],[450,550,650,750]):#range(nwins):I removed 750 from min,850 from max 
                wcnt=wcnt+1
                X[ii,wcnt,:,event_no]=np.mean(stc_cond.data[:,wcnt1:wcnt2],1)
#            X[ii,:,:,event_no]=np.transpose(stc_cond.data,[1,0]) 
#    X_list=[np.squeeze(x) for x in np.split(X, n_levels, axis=-1)]
#    factor_levels = [n_levels]  # number of levels in each factor
#    effects = 'A'  # this is the default signature for computing all effects
#    return_pvals = False          
#    def stat_fun(*args):  
#        return f_mway_rm(np.swapaxes(args, 1, 0), factor_levels=factor_levels,effects=effects, return_pvals=return_pvals)[0]

#            xx=np.asarray(args)
#            xx1=np.reshape(xx,(xx.shape[0],xx.shape[1],xx.shape[2]*xx.shape[3]))
#            xx2=np.transpose(xx1,[1,0,2])
#            return f_mway_rm(xx2, factor_levels=factor_levels,effects=effects, return_pvals=False)[0]

    source_space = grade_to_tris(5)
# as we only have one hemisphere we need only need half the connectivity
    print('Computing connectivity.')
    connectivity = spatial_tris_connectivity(source_space)

#    Now let's actually do the clustering. Please relax, on a small
#    notebook and one single thread only this will take a couple of minutes ...
    pthresh = p_threshold
#    max_step=1;
#    f_thresh = f_threshold_mway_rm(n_subjects, factor_levels, effects, pthresh)

#    To speed things up a bit we will ...
    n_permutations = 5000  # ... run fewer permutations (reduces sensitivity)
    fsave_vertices = [np.arange(10242), np.arange(10242)]
    if exclude_wbmedial:
        fname_label = label_path + '/' + 'toremove_wbspokes-lh.label'; labelL = mne.read_label(fname_label)
        fname_label = label_path + '/' + 'toremove_wbspokes-rh.label'; labelR = mne.read_label(fname_label)
        labelss=labelL+labelR
        bb=stc_cond.in_label(labelss)  
        nnl=np.in1d(fsave_vertices[0],bb.lh_vertno)
        nnr=np.in1d(fsave_vertices[1],bb.rh_vertno)
        spatial_exclude=np.hstack((fsave_vertices[0][nnl], fsave_vertices[0][nnr]+10242))
    if exclude_ROIs:
        fname_label='/imaging/rf02/Semnet/semnet4semloc//mask_labels_ATL_IFG_MTG_AG-lh.label'
        labelmask=mne.read_label(fname_label,subject='fsaverage')
        labelmask.values.fill(1.0)
        bb=stc_cond.in_label(labelmask)
        nnl=np.in1d(fsave_vertices[0],bb.lh_vertno)
        spatial_exclude=fsave_vertices[0][nnl].copy()

    print('Clustering.')
    t_threshold = -stats.distributions.t.ppf(p_threshold/2., n_subjects - 1)#dict(start=0, step=.1)#
    #t_threshold=2
    tail=0
    max_step=1
#    n_permutations=5000
    for refcond in np.asarray([0]):
#        refcond=thiscond+1
        X1=np.squeeze(X[:,:,:,thiscond]-X[:,:,:,refcond])
                            
        T_obs, clusters, cluster_p_values, H0 = clu = spatio_temporal_cluster_1samp_test(X1, connectivity=connectivity,  threshold=t_threshold, tail=tail, t_power=1, step_down_p=0.05, spatial_exclude=spatial_exclude, n_permutations=n_permutations, n_jobs=4, max_step=max_step)#, max_step=0)#spatial_exclude=exclude_ver, 
    
        #    T_obs, clusters, cluster_p_values, H0 = clu = \
        #        spatio_temporal_cluster_test(X_list, connectivity=connectivity, n_jobs=4,#max_step=max_step,
        #                                     threshold=f_thresh, stat_fun=stat_fun, spatial_exclude=spatial_exclude,
        #                                     n_permutations=n_permutations,
        #                                     buffer_size=None)
        #    Now select the clusters that are sig. at p < 0.05 (note that this value
        #    is multiple-comparisons corrected).
                                         
        p_thr=0.05
        good_cluster_inds = np.where(cluster_p_values <p_thr)[0]
        print cluster_p_values[good_cluster_inds]; print good_cluster_inds
        print cluster_p_values.min()
        ###############################################################################
        # Visualize the clusters
        
        print('Visualizing clusters.')
        
        vertices_to = [np.arange(10242), np.arange(10242)]
        #    tval_stc = mne.SourceEstimate(np.transpose(T_obs,[1,0]), vertices=fsave_vertices,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
        #    out_file3=uvttest_path + 'UVTtest_icomorphed_19subj_SemDec_pnt1_30ica_ConConds'
        #    tval_stc.save(out_file3)
        #    Now let's build a convenient representation of each cluster, where each
        #    cluster becomes a "time point" in the SourceEstimate
        if cluster_p_values.min()<=0.05:
            stc_all_cluster_vis = summarize_clusters_stc(clu, tstep=1e-3 * tstep1, vertices=fsave_vertices, subject='fsaverage', p_thresh=p_thr+0.0001)
            
            out_file1=out_path + 'ClusPer_abs_icomorphed_oldreg_clusterp'+str(p_threshold)[2:]+'_p'+str(p_thr)[2:]+'_18subj_LD_SL_1_48ica_'+event_names[thiscond]+'_'+event_names[refcond]+'_maxstep'+str(max_step)
            stc_all_cluster_vis.save(out_file1)
            
            Matx=np.zeros((20484,n_times))
            T=np.divide(T_obs,np.absolute(T_obs))
            for cc in range(good_cluster_inds.shape[0]):
              	Matx[clusters[good_cluster_inds[cc]][1],clusters[good_cluster_inds[cc]][0]]=T[clusters[good_cluster_inds[cc]][0],clusters[good_cluster_inds[cc]][1]]
            	
            
            #aa=Matx[clusters[good_cluster_inds[cc]][1],clusters[good_cluster_inds[cc]][0]]
            #bb=T_obs[clusters[good_cluster_inds[cc]][0],clusters[good_cluster_inds[cc]][1]]<0
            #aa[bb]=-1
            #Matx[clusters[good_cluster_inds[cc]][1],clusters[good_cluster_inds[cc]][0]]=aa
            
            
            matx_stc = mne.SourceEstimate(Matx, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
            out_file2=out_path + 'ClusPer_abs_sw_icomorphed_oldreg_clusterp'+str(p_threshold)[2:]+'_p'+str(p_thr)[2:]+'_18subj_LD_SL_1_48ica_'+event_names[thiscond]+'_'+event_names[refcond]+'_maxstep'+str(max_step)
            matx_stc.save(out_file2)
         
        
    	