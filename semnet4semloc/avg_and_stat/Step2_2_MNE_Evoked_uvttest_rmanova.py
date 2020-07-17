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
#sys.path.insert(1,'/imaging/rf02/Semnet/mne_python0.9')
#sys.path.insert(1,'/imaging/local/software/mne_python/latest')
sys.path.insert(1,'/imaging/local/software/mne_python/latest_v0.11')

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
from mne.stats import f_threshold_mway_rm, f_mway_rm, fdr_correction, spatio_temporal_cluster_test
from mne.datasets import sample
from mne.viz import mne_analyze_colormap
import os

###############################################################################
# Set parameters
out_path = '/imaging/rf02/Semnet/semnet4semloc/stc/permutation/evoked/' # root
uvttest_path = '/imaging/rf02/Semnet/semnet4semloc/stc/uvttest/evoked/rmanova/' # root

if not os.path.exists(out_path):
    os.makedirs(out_path)
if not os.path.exists(uvttest_path):
    os.makedirs(uvttest_path)
data_path = '/imaging/rf02/Semnet/'	# where subdirs for MEG data are
inv_path = '/imaging/rf02/Semnet/'
print (sys.argv)
p_inds = [0]
for ss in sys.argv[1:]:
   p_inds.append( int( ss ) )
label_path = '/imaging/rf02/TypLexMEG/fsaverage/label'

#subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19] # removed
#p_inds=[0,1,2,3,4,5,6,7,8,9,10]
p_list=[0.00005,0.05,0.045,0.04,0.03,0.025,0.02,0.01,0.008,0.005,0.002,0.001, 0.0008,0.0005,0.0002,0.0001,0.00005,0.00001]#0.001,0.0009,0.0008,0.0007,0.0006,0.0005]#
#print ("subject_inds:")
#print (subject_inds)
print ("No rejection")
b='gamma'
list_all =  ['/meg16_0032/160218/', #1  '/meg16_0030/160216/', #0
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

print ("ll:")
print (ll)

n_subjects=18
# so we have X: observations * time * vertices in each condition
#X_list=range(2)
semtasks=['SemDec','LD']#
event_names = ['Emotional', 'Concrete']#, ,'Pwordc']'Neutral', 'Emotional',
effect_names=['task','contrast']
all_effects=['A','B']

n_levels=len(semtasks)
n_factors=len(event_names)
n_times=len(list(np.arange(350,751,100)))
tmin1=50
tstep1=100
X=np.zeros((n_subjects,n_levels*n_factors,20484,n_times))#np.zeros((n_subjects,n_times,20484,n_levels))
nwins=5
#Xmean=np.zeros((n_subjects,nwins,20484,n_levels))
for p_threshold in ll:
    
    for ii, meg in enumerate(list_all):
        print (ii)
        tecnt=-1
        for tcnt,task_name in enumerate(semtasks):#task is A
            print(b)
            for evcnt, event_name in enumerate(event_names):#contrast is B
                tecnt=tecnt+1
                fname_in = data_path + meg + 'firstMorphed_ico_oldreg_'+task_name+'_SL_1_48ica_'+event_name+'_Source_Evoked_m300_600'#_'+b#'firstMorphed_ico_SemDec_ica_'+event_names[event_no]+'_Source_Evoked_mnrsmpl_50_550_100ms_0ov'#'_Source_Evoked_m300_600'#'_Source_Evoked_mnrsmpl_50_550_100ms_0ov'
                stc_cond = mne.read_source_estimate(fname_in)
    #            stc_cond.resample(100)
    #            stc_cond.crop(0.050,0.450)
                wcnt=-1
                for wcnt1,wcnt2 in zip(list(np.arange(350,751,100)),list(np.arange(450,851,100))):#range(nwins):
                    print wcnt1,wcnt2
                    wcnt=wcnt+1
                    X[ii,tecnt,:,wcnt]=np.mean(stc_cond.data[:,wcnt1:wcnt2],1)

                #X[ii,:,:,event_no]=np.transpose(stc_cond.data,[1,0]) #[:,350:650]
    X1=X.copy()#[:,:,:,np.array([0,2,1,3])].copy()
    
    #=[np.squeeze(x) for x in np.split(X1, n_levels*n_factors, axis=-1)]
    factor_levels = [n_levels,n_levels]  # number of levels in each factor
    for this_effect,effect_name in zip(all_effects,effect_names)
        #effects = 'B'  # A*B is the default signature for computing all effects, A here is task effect, B contrast
        return_pvals = True         
        #def stat_fun(*args):  
        #    return f_mway_rm(np.swapaxes(args, 1, 0), factor_levels=factor_levels,effects=effects, return_pvals=return_pvals)[0]
        #thisX=np.swapaxes(X1,1,3)#subxcondsxtime
        f=np.zeros((X1.shape[2],X1.shape[3]))
        for xii in range(X1.shape[3]):#over time
            f[:,xii], p= f_mway_rm(X1[:,:,:,xii], factor_levels=factor_levels,effects=this_effect, return_pvals=return_pvals)
        matx_stc = mne.SourceEstimate(f, vertices=vertices_avg,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
        out_file2=uvttest_path + 'UVttest_rmANOVA_Evoked_icomorphed_oldreg_18subj_SDvsLD_'+effect_name
        matx_stc.save(out_file2)
    
#    ###############################################################################
#    # Visualize the clusters
#    
#    print('Visualizing clusters.')
#    
#    #    tval_stc = mne.SourceEstimate(np.transpose(T_obs,[1,0]), vertices=fsave_vertices,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
#    #    out_file3=uvttest_path + 'UVTtest_icomorphed_19subj_SemDec_pnt1_30ica_ConConds'
#    #    tval_stc.save(out_file3)
#    #    Now let's build a convenient representation of each cluster, where each
#    #    cluster becomes a "time point" in the SourceEstimate
#    if len(cluster_p_values)>0:
#        if cluster_p_values.min()<0.05:
#            stc_all_cluster_vis = summarize_clusters_stc(clu, tstep=1e-3 * tstep1, vertices=fsave_vertices, subject='fsaverage', p_thresh=p_thr+0.0001)
#            
#            out_file1=out_path + 'ClusPer_sEvoked_icomorphed_newreg_clusterp'+str(p_threshold)[2:]+'_p'+str(p_thr)[2:]+'_19subj_SemDec_pnt1_30ica_'+event_names[0]+'_'+event_names[1]+'gammamu_maxstep'+str(max_step)
#            stc_all_cluster_vis.save(out_file1)
#            
#            Matx=np.zeros((20484,n_times))
#            T=np.divide(T_obs,np.absolute(T_obs))
#            for cc in range(good_cluster_inds.shape[0]):
#              	Matx[clusters[good_cluster_inds[cc]][1],clusters[good_cluster_inds[cc]][0]]=T[clusters[good_cluster_inds[cc]][0],clusters[good_cluster_inds[cc]][1]]
#            	
#            
#            #aa=Matx[clusters[good_cluster_inds[cc]][1],clusters[good_cluster_inds[cc]][0]]
#            #bb=T_obs[clusters[good_cluster_inds[cc]][0],clusters[good_cluster_inds[cc]][1]]<0
#            #aa[bb]=-1
#            #Matx[clusters[good_cluster_inds[cc]][1],clusters[good_cluster_inds[cc]][0]]=aa
#            
#            
#            matx_stc = mne.SourceEstimate(Matx, vertices=vertices_avg,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
#            out_file2=out_path + 'ClusPer_sEvoked_sw_icomorphed_newreg_clusterp'+str(p_threshold)[2:]+'_p'+str(p_thr)[2:]+'_19subj_SemDec_pnt1_30ica_'+event_names[0]+'_'+event_names[1]+'gammamu_maxstep'+str(max_step)
#            matx_stc.save(out_file2)
#        
        
     
        
    	
