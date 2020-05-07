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
subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18] # removed
nsubs=len(subject_inds)
out_path = '/imaging/rf02/Semnet/stc/permutation/evoked/bands/newfilt/newregsigned/taugamma_newev/nsubs_'+str(nsubs)+'/' # root
uvttest_path = '/imaging/rf02/Semnet/stc/uvttest/evoked/bands/newfilt/newregsigned/taugamma_newev/nsubs_'+str(nsubs)+'/' # root

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


#p_inds=[0,1,2,3,4,5,6,7,8,9,10]
p_list=[0.05, 0.025, 0.01, 0.008, 0.005]#0.001,0.0009,0.0008,0.0007,0.0006,0.0005]#
print "subject_inds:"
print subject_inds
print "No rejection"

list_all =  ['/meg16_0030/160216/', #0
            '/meg16_0032/160218/', #1
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
            '/meg16_0122/160707/', #22 LD
#            '/meg16_0123/160708/', #23 LD
            '/meg16_0125/160712/', #24 LD
            ]

# subjects names used for MRI data
subjects=  ['MRI_meg16_0030' ,#0
            'MRI_meg16_0032' ,#1
            'MRI_meg16_0034' ,#3
            'MRI_meg16_0035' ,#4
            'MRI_meg16_0042' ,#7
            'MRI_meg16_0045' ,#8
            'MRI_meg16_0052' ,#10
            'MRI_meg16_0056' ,#11
    	       'MRI_meg16_0069' ,#12
 	       'MRI_meg16_0070' ,#13
            'MRI_meg16_0072' ,#15
            'MRI_meg16_0073' ,#16
            'MRI_meg16_0075' ,#17
            'MRI_meg16_0078' ,#18
            'MRI_meg16_0082' ,#19
            'MRI_meg16_0086' ,#20
            'MRI_meg16_0097' ,#21
            'MRI_meg16_0122' ,#22
#            'MRI_meg16_0123' ,#23
            'MRI_meg16_0125' ,#24
            ]

ll = []
for ss in p_inds:
 ll.append(p_list[ss])

print "ll:"
print ll

n_subjects=len(subject_inds)
#stim_delay = 0.034 # delay in s ((NOTE! not in the example, taken from Olaf's script. Rezvan))
tmin, tmax = -0.3, 0.6
# so we have X: observations * time * vertices in each condition
#X_list=range(2)
event_names = ['Visual','Hear']# 'Visual']#'Hear']#, 'Neutral', 'Emotional','Pwordc']
bands=['alpha','gamma']#['delta','theta','alpha','beta','gamma']
#thiscond=3
#refcond=1
n_levels=len(event_names)
n_factors=len(bands)
n_times=len(list(np.arange(350,751,25)))
tmin1=50.
tstep1=25.
#nwins=5
X=np.zeros((n_subjects,n_times,20484,n_levels*n_factors))
for p_threshold in ll:
#    evbncnt=-1
    for ii, meg in enumerate(list_all):
        print subject_inds[ii]
        for bcnt,b in enumerate(bands):
            print b
            for evcnt, event_name in enumerate(event_names):                
                fname_in = data_path + meg + 'firstMorphed_SemDec_filtnew_pnt1_48ica_'+event_name+'_Source_signedEvoked_m300_600_'+b#'firstMorphed_ico_SemDec_pnt1_30ica_'+event_names[event_no]+'_Source_Evoked_mnrsmpl_50_550_100ms_0ov'#'_Source_Evoked_m300_600'#'_Source_Evoked_mnrsmpl_50_550_100ms_0ov'

                stc_cond = read_source_estimate(fname_in)
        
                method = "MNE"  # use dSPM method (could also be MNE or sLORETA)                
                snr = 3.0#np.mean(ediv)
                print snr
                lambda2 = 1.0 / snr ** 2
                subject_from = subjects[subject_inds[ii]]
                subject_to = 'fsaverage'
                vertices_avg = [np.arange(10242), np.arange(10242)]      
    
                #    stc_list=range(len(epochs_list))
                srcin='/imaging/rf02/TypLexMEG/fsaverage/bem/fsaverage_dist-ico-5-src.fif'
                src_avg = mne.read_source_spaces(srcin)
                wcnt=-1
            for wcnt1,wcnt2 in zip(list(np.arange(350,751,25)),list(np.arange(360,761,25))):#range(nwins):
#                print wcnt1,wcnt2
                wcnt=wcnt+1
                X[ii,wcnt,:,bcnt+evcnt]=np.mean(stc_cond.data[:,wcnt1:wcnt2],1)

#            stc_cond.resample(60)
#            stc_cond.crop(0.050,0.450)
#            X[ii,:,:,event_no]=np.transpose(stc_cond.data,[1,0]) 
            
#            X[ii,:,:,event_no]=np.transpose(stc_cond.data,[1,0])
    X1=X.copy()#[:,:,:,np.array([0,2,1,3])].copy()
    X_list=[np.squeeze(x) for x in np.split(X1, n_levels*n_factors, axis=-1)]
    factor_levels = [n_factors,n_factors]  # number of levels in each factor
    effects = 'A:B'  # this is the default signature for computing all effects
    return_pvals = False          
    def stat_fun(*args):  
        return f_mway_rm(np.swapaxes(args, 1, 0), factor_levels=factor_levels,effects=effects, return_pvals=return_pvals)[0]

    source_space = grade_to_tris(5)
# as we only have one hemisphere we need only need half the connectivity
    print('Computing connectivity.')
    connectivity = spatial_tris_connectivity(source_space)

#    Now let's actually do the clustering. Please relax, on a small
#    notebook and one single thread only this will take a couple of minutes ...
    pthresh = p_threshold
#    max_step=1;
    f_thresh = f_threshold_mway_rm(n_subjects, factor_levels, effects, pthresh)

#    To speed things up a bit we will ...
    n_permutations = 10000  # ... run fewer permutations (reduces sensitivity)
    fname_label = label_path + '/' + 'toremove_wbspokes-lh.label'; labelL = mne.read_label(fname_label)
    fname_label = label_path + '/' + 'toremove_wbspokes-rh.label'; labelR = mne.read_label(fname_label)
    labelss=labelL+labelR
    bb=stc_cond.in_label(labelss)
    fsave_vertices = [np.arange(10242), np.arange(10242)]
    nnl=np.in1d(fsave_vertices[0],bb.lh_vertno)
    nnr=np.in1d(fsave_vertices[1],bb.rh_vertno)
    spatial_exclude=np.hstack((fsave_vertices[0][nnl], fsave_vertices[0][nnr]+10242))
    print('Clustering.')
#    t_threshold = -stats.distributions.t.ppf(p_threshold/2., n_subjects - 1)#dict(start=0, step=.1)#
    #t_threshold=2
    tail=0
    max_step=1;
#    n_permutations=5000
#    for thiscond, refcond in zip(list(np.array([0, 0, 1])),list(np.array([1, 2, 2]))):
##        print thiscond, refcond
##        refcond=thiscond+1
##        if thiscond<3:
##            X1=np.squeeze(X[:,:,:,thiscond]-X[:,:,:,refcond])
##        else:
##            X1=np.squeeze(X[:,:,:,refcond]-X[:,:,:,thiscond])
#        X1=np.squeeze(X[:,:,:,thiscond]-X[:,:,:,refcond])   
#        ##uvttest
#        Matx_uvttest=np.transpose(X1.copy(),[0,2,1])#np.zeros((n_subjects,20484,n_times))
#        MT=np.zeros((20484,n_times))
#    
#        for cnt in range(Matx_uvttest.shape[2]):
#            print cnt
#            MT[:,cnt]= mne.stats.ttest_1samp_no_p(Matx_uvttest[:,:,cnt], sigma=1e-3, method='relative')
#        #MT[np.abs(MT)<t_threshold]=0
#    
#        tval_stc = mne.SourceEstimate(MT, vertices=vertices_avg,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
#        out_file3=uvttest_path + 'UVTtest_icomorphed_newreg_19subj_SemDec_m300_600_100ms_pnt1_30ica_'+event_names[thiscond]+'_'+event_names[refcond]#+'_'+b
#        tval_stc.save(out_file3)
#                
#        T_obs, clusters, cluster_p_values, H0 = clu = spatio_temporal_cluster_1samp_test(X1, connectivity=connectivity,  threshold=t_threshold, tail=tail, t_power=1, step_down_p=0.05, spatial_exclude=spatial_exclude, n_permutations=n_permutations, n_jobs=4)#, max_step=0)#, max_step=0)#spatial_exclude=exclude_ver, 
    
    T_obs, clusters, cluster_p_values, H0 = clu = \
        spatio_temporal_cluster_test(X_list, connectivity=connectivity, n_jobs=16,#step_down_p=0.05,#max_step=max_step,
                                     threshold=f_thresh, stat_fun=stat_fun, spatial_exclude=spatial_exclude,
                                     n_permutations=n_permutations,
                                     buffer_size=None)
#        Now select the clusters that are sig. at p < 0.05 (note that this value
#        is multiple-comparisons corrected).
                                     
    p_thr=0.05
    good_cluster_inds = np.where(cluster_p_values <p_thr)[0]
    print cluster_p_values[good_cluster_inds]; print good_cluster_inds
#        print cluster_p_values.min()
    ###############################################################################
    # Visualize the clusters
    
    print('Visualizing clusters.')
    
    #    tval_stc = mne.SourceEstimate(np.transpose(T_obs,[1,0]), vertices=fsave_vertices,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
    #    out_file3=uvttest_path + 'UVTtest_icomorphed_19subj_SemDec_pnt1_30ica_ConConds'
    #    tval_stc.save(out_file3)
    #    Now let's build a convenient representation of each cluster, where each
    #    cluster becomes a "time point" in the SourceEstimate
    if len(cluster_p_values)>0:
        if cluster_p_values.min()<0.05:
            stc_all_cluster_vis = summarize_clusters_stc(clu, tstep=1e-3 * tstep1, vertices=fsave_vertices, subject='fsaverage', p_thresh=p_thr+0.0001)
            
            out_file1=out_path + 'ClusPer_sEvoked_icomorphed_newreg_clusterp'+str(p_threshold)[2:]+'_p'+str(p_thr)[2:]+'_20subj_SemDec_pnt1_30ica_'+event_names[0]+'_'+event_names[1]+'gammatau_maxstep'+str(max_step)
            stc_all_cluster_vis.save(out_file1)
            
            Matx=np.zeros((20484,n_times))
            T=np.divide(T_obs,np.absolute(T_obs))
            for cc in range(good_cluster_inds.shape[0]):
              	Matx[clusters[good_cluster_inds[cc]][1],clusters[good_cluster_inds[cc]][0]]=T[clusters[good_cluster_inds[cc]][0],clusters[good_cluster_inds[cc]][1]]
            	
            
            #aa=Matx[clusters[good_cluster_inds[cc]][1],clusters[good_cluster_inds[cc]][0]]
            #bb=T_obs[clusters[good_cluster_inds[cc]][0],clusters[good_cluster_inds[cc]][1]]<0
            #aa[bb]=-1
            #Matx[clusters[good_cluster_inds[cc]][1],clusters[good_cluster_inds[cc]][0]]=aa
            
            
            matx_stc = mne.SourceEstimate(Matx, vertices=vertices_avg,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
            out_file2=out_path + 'ClusPer_sEvoked_sw_icomorphed_newreg_clusterp'+str(p_threshold)[2:]+'_p'+str(p_thr)[2:]+'_20subj_SemDec_pnt1_30ica_'+event_names[0]+'_'+event_names[1]+'gammatau_maxstep'+str(max_step)
            matx_stc.save(out_file2)
        
        
     
        
    	