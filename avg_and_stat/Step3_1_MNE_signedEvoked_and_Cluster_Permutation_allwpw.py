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
out_path = '/imaging/rf02/Semnet/stc/permutation/evoked/pnt1_30/newregsigned/allwpw/' # root
uvttest_path = '/imaging/rf02/Semnet/stc/uvttest/evoked/pnt1_30/newregsigned/allwpw/' # root
avg_path= '/imaging/rf02/Semnet/stc/GrandAverage/evoked/pnt1_30/newregsigned/allwpw/'

if not os.path.exists(out_path):
    os.makedirs(out_path)
if not os.path.exists(uvttest_path):
    os.makedirs(uvttest_path)
if not os.path.exists(avg_path):
    os.makedirs(avg_path)
data_path = '/imaging/rf02/Semnet/'	# where subdirs for MEG data are
inv_path = '/imaging/rf02/Semnet/'
print sys.argv
p_inds = []
for ss in sys.argv[1:]:
   p_inds.append( int( ss ) )
label_path = '/imaging/rf02/TypLexMEG/fsaverage/label'

subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19] # removed
#p_inds=[0,1,2,3,4,5,6,7,8,9,10]
p_list=[0.05,0.025,0.01,0.008,0.005]#,0.002,0.001, 0.0008,0.0005,0.0002,0.0001,0.00005,0.00001]#0.001,0.0009,0.0008,0.0007,0.0006,0.0005]#
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
            '/meg16_0122/160707/', #22 
            '/meg16_0125/160712/', #24 
            ]



# subjects names used for MRI data
subjects=  ['MRI_meg16_0030' ,#0
            'MRI_meg16_0032' ,#1
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

n_subjects=19
stim_delay = 0.034 # delay in s ((NOTE! not in the example, taken from Olaf's script. Rezvan))
tmin, tmax = -0.3, 0.6
# so we have X: observations * time * vertices in each condition
#X_list=range(2)
event_names = ['Conwords','Pwordc']# 'Visual']#'Hear']#, 'Neutral', 'Emotional','Pwordc']
#thiscond=3
refcond=1
n_levels=len(event_names)
n_times=len(list(np.arange(350,751,25)))
tmin1=0
tstep1=25.
#nwins=5
X=np.zeros((n_subjects,n_times,20484,n_levels))
for p_threshold in ll:
    
    for ii, meg in enumerate(list_all):
        print subject_inds[ii]
        #cov_path=data_path+meg+'noise_cov'
        #noise_cov=mne.read_cov(cov_path)
        raw_fname = data_path + meg+ 'SemDec_blocks_tsss_filt_pnt1_30_ica_raw.fif' #clean_ssp_
        event_fname = data_path + meg + 'SemDec_blocks_tsss_filt_pnt1_30_ica_raw-eve.fif'
        subjects_dir = data_path 
        
        tmin = -0.3
        tmax = 0.6  # Use a lower tmax to reduce multiple comparisons
        #   Setup for reading the raw data
        raw = io.Raw(raw_fname)
        events = mne.read_events(event_fname)
        stim_delay = 0.034 # delay in s
        events[:,0] = events[:,0]+np.round( raw.info['sfreq']*stim_delay )
        ###############################################################################
        # Read epochs for all channels, removing a bad one
        print "epochs"
        picks = mne.pick_types(raw.info, eeg=True, meg=True, eog=False, exclude='bads')
        print "Epoching"
        if subject_inds[ii] in np.array([3,5]):#np.array([4,8,9]):
            reject = dict(eeg=120e-6, grad=200e-12, mag=4e-12)
        else:
            reject = dict(eeg=100e-6, grad=200e-12, mag=4e-12)#, eog=150e-6)
        event_ids_w = {'visual': 1, 'hear': 2, 'hand': 3}#, 'neutral': 4, 'emotional': 5,'pwordc': 6}
        epochs_w = mne.Epochs(raw, events, event_ids_w, tmin, tmax, picks=picks, proj=True, baseline=(None, 0), reject=reject)
        epochs_w.drop_bad_epochs()
        event_ids_pw = {'pwordc': 6}
        epochs_pw = mne.Epochs(raw, events, event_ids_pw, tmin, tmax, picks=picks, proj=True, baseline=(None, 0), reject=reject)
        epochs_pw.drop_bad_epochs()
        #    for eegcnt in range(71):
        #        if eegcnt<10:
        #            thiseegbad=sum(epochs.drop_log,[]).count('EEG00'+str(eegcnt))
        #            if thiseegbad>=90:
        #                raw.info['bads'].append(u'EEG00'+str(eegcnt))                
        #        else:
        #            thiseegbad=sum(epochs.drop_log,[]).count('EEG0'+str(eegcnt))
        #            if thiseegbad>=90:
        #                raw.info['bads'].append(u'EEG0'+str(eegcnt))
        picks = mne.pick_types(raw.info, eeg=True, meg=True, eog=False, exclude='bads')
        epochs_list=[epochs_w,epochs_pw]#range(n_levels)
#        for eid,event_id in enumerate(list(np.array([1,2,3,6]))):
##                event_id = eid+1  
#            epochs_list[eid]=epochs[epochs.events[:,2]==event_id]
#        equalize_epoch_counts(epochs_list,method='mintime')

               
               
               
    	#    Equalize trial counts to eliminate bias (which would otherwise be
    	#    introduced by the abs() performed below)
    
    	###############################################################################
        print "Transform to source space"
        
        fname_inv = data_path + meg + 'InvOp_ico5newreg_fft_pnt1_30_clean_ica_EMEG-inv.fif'
        
        method = "MNE"  # use dSPM method (could also be MNE or sLORETA)
        inverse_operator = read_inverse_operator(fname_inv)
        
        snr = 3.0#np.mean(ediv)
        print snr
        lambda2 = 1.0 / snr ** 2
        subject_from = subjects[subject_inds[ii]]
        subject_to = 'fsaverage'
        vertices_avg = [np.arange(10242), np.arange(10242)]      
    
        #    stc_list=range(len(epochs_list))
        srcin='/imaging/rf02/TypLexMEG/fsaverage/bem/fsaverage_dist-ico-5-src.fif'
        src_avg = mne.read_source_spaces(srcin)
        tminl=0
        tstepl=1
        thislbl_data=np.ones((len(vertices_avg[0])+len(vertices_avg[1]),1))
        thisstc_tolabel = mne.SourceEstimate(thislbl_data, vertices=vertices_avg,tmin=1e-3 * tminl, tstep=1e-3 * tstepl, subject='fsaverage')
        thisstc_label=mne.stc_to_label(thisstc_tolabel, src_avg)
        thisscale=-np.hstack((mne.label.label_sign_flip(thisstc_label[0],src_avg),mne.label.label_sign_flip(thisstc_label[1],src_avg)))
        for evcnt in range(len(epochs_list)):
            this_evoked = epochs_list[evcnt].copy().average()
            this_condition = apply_inverse(this_evoked, inverse_operator, lambda2, method, pick_ori="normal")
            this_morphmat=mne.compute_morph_matrix(subject_from, subject_to, this_condition.vertices, vertices_avg, subjects_dir=data_path)
            this_stc=this_condition.morph_precomputed(subject_to,vertices_avg,this_morphmat, subject_from)
            tstepds=1/raw.info['sfreq']
            predata1=this_stc.data.copy()
            ntimes=len(this_stc.times)
            #    label_data1=np.tile(thisscale,[ntimes,1]).transpose()*label_predata
            signeddata=predata1.copy()#np.tile(thisscale[:,np.newaxis],(1,predata1.shape[1]))*predata1
            stc_signeddata = mne.SourceEstimate(signeddata, vertices=vertices_avg,tmin=1e-3 * tmin, tstep=tstepds, subject='fsaverage')
            wcnt=-1
            for wcnt1,wcnt2 in zip(list(np.arange(350,751,25)),list(np.arange(375,776,25))):#range(nwins):
                wcnt=wcnt+1
                X[ii,wcnt,:,evcnt]=np.mean(stc_signeddata.data[:,wcnt1:wcnt2],1)

#            stc_cond.resample(60)
#            stc_cond.crop(0.050,0.450)
#            X[ii,:,:,event_no]=np.transpose(stc_cond.data,[1,0]) 
            
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
    fname_label = label_path + '/' + 'toremove_wbspokes-lh.label'; labelL = mne.read_label(fname_label)
    fname_label = label_path + '/' + 'toremove_wbspokes-rh.label'; labelR = mne.read_label(fname_label)
    labelss=labelL+labelR
    bb=stc_signeddata.in_label(labelss)
    fsave_vertices = [np.arange(10242), np.arange(10242)]
    nnl=np.in1d(fsave_vertices[0],bb.lh_vertno)
    nnr=np.in1d(fsave_vertices[1],bb.rh_vertno)
    spatial_exclude=np.hstack((fsave_vertices[0][nnl], fsave_vertices[0][nnr]+10242))
    print('Clustering.')
    t_threshold = -stats.distributions.t.ppf(p_threshold/2., n_subjects - 1)#dict(start=0, step=.1)#
    #t_threshold=2
    tail=0
    max_step=5;
#    n_permutations=5000
    for condavg in range(n_levels):
        Xavg=np.squeeze(np.mean(X[:,:,:,condavg],0)).T
        avg_stc = mne.SourceEstimate(Xavg, vertices=vertices_avg,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
        out_file4=avg_path + 'Grandavg_icomorphed_newreg_19subj_SemDec_m300_600_25ms_pnt1_30ica_'+event_names[condavg]#+'_'+b
        avg_stc.save(out_file4)
        
        
    for thiscond in range(1):
#        refcond=thiscond+1
        if thiscond<3:
            X1=np.squeeze(X[:,:,:,thiscond]-X[:,:,:,refcond])
        else:
            X1=np.squeeze(X[:,:,:,refcond]-X[:,:,:,thiscond])
            
        ##uvttest
        Matx_uvttest=np.transpose(X1.copy(),[0,2,1])#np.zeros((n_subjects,20484,n_times))
        MT=np.zeros((20484,n_times))
    
        for cnt in range(Matx_uvttest.shape[2]):
            print cnt
            MT[:,cnt]= mne.stats.ttest_1samp_no_p(Matx_uvttest[:,:,cnt], sigma=1e-3, method='relative')
        #MT[np.abs(MT)<t_threshold]=0
    
        tval_stc = mne.SourceEstimate(MT, vertices=vertices_avg,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
        out_file3=uvttest_path + 'UVTtest_icomorphed_newreg_19subj_SemDec_m300_600_25ms_pnt1_30ica_'+event_names[thiscond]+'_'+event_names[refcond]#+'_'+b
        tval_stc.save(out_file3)
                
        T_obs, clusters, cluster_p_values, H0 = clu = spatio_temporal_cluster_1samp_test(X1, connectivity=connectivity,  threshold=t_threshold, tail=tail, t_power=1, step_down_p=0.05, spatial_exclude=spatial_exclude, n_permutations=n_permutations, n_jobs=8)#, max_step=0)#, max_step=0)#spatial_exclude=exclude_ver, 
    
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
                
                out_file1=out_path + 'ClusPer_abs_icomorphed_newreg_clusterp'+str(p_threshold)[2:]+'_p'+str(p_thr)[2:]+'_19subj_SemDec_pnt1_30ica_'+event_names[thiscond]+'_'+event_names[refcond]+'_maxstep'+str(max_step)
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
                out_file2=out_path + 'ClusPer_abs_sw_icomorphed_newreg_clusterp'+str(p_threshold)[2:]+'_p'+str(p_thr)[2:]+'_19subj_SemDec_pnt1_30ica_'+event_names[thiscond]+'_'+event_names[refcond]+'_maxstep'+str(max_step)
                matx_stc.save(out_file2)
                
                
         
#        
#    	