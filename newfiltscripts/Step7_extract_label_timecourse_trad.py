# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 14:55:11 2015

@author: rf02
"""

# Author: Martin Luessi <mluessi@nmr.mgh.harvard.edu>
#
# License: BSD (3-clause)

print(__doc__)
import sys

#sys.path.append('/imaging/local/software/python_packages/nibabel/1.3.0')
#sys.path.append('/imaging/local/software/python_packages/pysurfer/v0.4')
#sys.path.append('/imaging/local/software/python_packages/scikit-learn/v0.14.1')
#sys.path.append('/imaging/local/software/mne_python/latest')

# for qsub
# (add ! here if needed)
sys.path.append('/imaging/local/software/EPD/latest/x86_64/bin/python')
sys.path.insert(1,'/imaging/rf02/TypLexMEG/mne_python0.9')
#sys.path.insert(1,'/imaging/local/software/mne_python/latest_v0.9')
sys.path.append('/imaging/local/software/python_packages/scikit-learn/v0.16.1')
sys.path.insert(1,'/imaging/local/software/python_packages/scikit-learn/v0.13.1/lib.linux-x86_64-2.7/sklearn/externals')
import joblib
###


import os
import numpy as np
import mne
import scipy
from scipy.stats.mstats import zscore
#import sklearn
from mne.io import Raw
from mne.epochs import equalize_epoch_counts
from mne.minimum_norm import (read_inverse_operator, psf_ctf_28Oct15, apply_inverse)
from mne.stats import f_threshold_mway_rm, f_mway_rm
from mne import (io, spatial_tris_connectivity, grade_to_tris)
from mne.epochs import equalize_epoch_counts
from mne.stats import (permutation_cluster_test,summarize_clusters_stc)
from mne.minimum_norm import apply_inverse, read_inverse_operator
from mne.stats import f_threshold_mway_rm, f_mway_rm, fdr_correction, spatio_temporal_cluster_test
from mne.datasets import sample
from mne.viz import mne_analyze_colormap
import matplotlib.pyplot as plt
import copy
main_path = '/imaging/rf02/Semnet/'
data_path = '/imaging/rf02/Semnet/'	# where subdirs for MEG data are
event_path = '/imaging/rf02/Semnet'
os.chdir(data_path)
label_path = '/imaging/rf02/Semnet/fsaverage/label/'
inv_path = '/imaging/rf02/Semnet/'
subjects_dir = data_path
print sys.argv
subject_inds = []
for ss in sys.argv[1:]:
   subject_inds.append( int( ss ) )

subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]
print "subject_inds:"
print subject_inds
print "No rejection"

list_all =  ['/meg16_0030/160216/', #0
            '/meg16_0032/160218/', #1
#            '/meg16_0033/160218/', #2
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
            '/meg16_0123/160708/', #23 LD
            '/meg16_0125/160712/', #24 LD
            ]

# subjects names used for MRI data
subjects=  ['MRI_meg16_0030' ,#0
            'MRI_meg16_0032' ,#1
            #            'MRI_meg16_0033' ,#2
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
            'MRI_meg16_0123' ,#23
            'MRI_meg16_0125' ,#24
            ]

ll = []
for ss in subject_inds:
 ll.append(list_all[ss])

print "ll:"
print ll

labellist_path = [label_path+'vis_final-lh.label',label_path+'auditory_final_manual3-lh.label',label_path+'hand_c2-lh.label',label_path+'ATL_manual-lh.label',label_path+'AG_manual3-lh.label',label_path+'lh.parsorbitalis.label']#[ 'lh.posterior_fusiform','rh.posterior_fusiform', 'lh.lateraloccipital','rh.lateraloccipital','lh.precentral','rh.precentral', 'lh.ATL_latmed','rh.ATL_latmed']#labellist = [ 'leftATL-lh','rightATL-rh']#, 'medialtempright-rh','medialtempleft-lh']'lh.postcentral','rh.postcentral',
labellist1=[mne.read_label(label) for label in labellist_path[:3]]
for thisl in labellist1:
    thisl.values.fill(1.0)
fsaverage_path='/imaging/rf02/Semnet'
labelss = mne.read_labels_from_annot(subject='fsaverage', parc='rez_aparc_19oct16', subjects_dir=fsaverage_path)
#labelss.remove(labelss[-3])
#labelss.remove(labelss[-2])
labelss.remove(labelss[-1])
for lbl in labelss:
    lbl.values.fill(1.0)
label_names=[label.name for label in labelss]
label_index=np.array([label_names.index('middle_temporal-lh')])#,label_names.index('AG_hsmetafinal5-lh')])

labellist_hubs=[mne.read_label(label) for label in labellist_path[3:]]+[labelss[ii] for ii in label_index]#[labelss[ii] for ii in label_index]
for thisl in labellist_hubs:
    thisl.values.fill(1.0)
    
    
###RH
#labellist_path = [label_path+'vis_primary2-rh.label',label_path+'auditory_primary1-rh.label',label_path+'hand_c-rh.label',label_path+'ATL_manual-rh.label',label_path+'AG_manual-rh.label',label_path+'rh.parsorbitalis.label']#[ 'lh.posterior_fusiform','rh.posterior_fusiform', 'lh.lateraloccipital','rh.lateraloccipital','lh.precentral','rh.precentral', 'lh.ATL_latmed','rh.ATL_latmed']#labellist = [ 'leftATL-lh','rightATL-rh']#, 'medialtempright-rh','medialtempleft-lh']'lh.postcentral','rh.postcentral',
#labellist1=[mne.read_label(label) for label in labellist_path[:3]]
#for thisl in labellist1:
#    thisl.values.fill(1.0)
#fsaverage_path='/imaging/rf02/Semnet'
#labelss = mne.read_labels_from_annot(subject='fsaverage', parc='rez_aparc_19oct16', subjects_dir=fsaverage_path)
##labelss.remove(labelss[-3])
##labelss.remove(labelss[-2])
#labelss.remove(labelss[-1])
#for lbl in labelss:
#    lbl.values.fill(1.0)
#label_names=[label.name for label in labelss]
#label_index=np.array([label_names.index('middle_temporal-rh')])#,label_names.index('AG_hsmetafinal5-lh')])
#
#labellist_hubs=[mne.read_label(label) for label in labellist_path[3:]]+[labelss[ii] for ii in label_index]#[labelss[ii] for ii in label_index]
#for thisl in labellist_hubs:
#    thisl.values.fill(1.0)
    
##############
#print labellist
#mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/lh.supramarginal.label', supramarginal)
#
#posterior_inferiorfrontal=labelss[label_name.index('posterior_inferiorfrontal-lh')]
#mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/lh.posterior_inferiorfrontal.label', posterior_inferiorfrontal)
#
#middle_temporal=labelss[label_name.index('middle_temporal-lh')]
#mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/lh.middle_temporal.label', middle_temporal)
event_names = ['Visual', 'Hear', 'Hand']#,'Pwordc']#, 'Neutral', 'Emotional','Pwordc']
sfrq=5.
n_subjects=len(list_all)
n_times=int(901/sfrq)
ntimes=901
which_conds=np.array([0,1,2])
labellist_spokes=[labellist1[wc] for wc in which_conds]
labellist=labellist_hubs+labellist_spokes
n_levels_facts=len(which_conds)*len(labellist)
#X=np.zeros((n_subjects,n_times,n_levels_facts))
nconds=len(which_conds)
nrois=len(labellist)
tcmat_normf=np.zeros((n_subjects,nrois,ntimes,nconds))
tcmat_nonm=np.zeros((n_subjects,nrois,ntimes,nconds))
nevents=len(event_names)
print labellist
out_path=data_path+'ROI_analysis/'
if not os.path.exists(out_path):
    os.makedirs(out_path)


#event_no=0
ii=-1
for index, meg in enumerate(ll):
    print subject_inds[index]
    subject_no=subject_inds[index]
    print "compute SNR"
    raw_fname = data_path + meg+ 'SemDec_blocks_tsss_filtnew_pnt1_30_ica_raw.fif' #clean_ssp_
    event_fname = data_path + meg + 'SemDec_blocks_tsss_filtnew_pnt1_30_ica_raw-eve.fif'
    subjects_dir = data_path 
    
    tmin = -0.3
    sonset=int(abs(tmin*1000)/sfrq)
    tmax = 0.6  # Use a lower tmax to reduce multiple comparisons
    #   Setup for reading the raw data
    raw = Raw(raw_fname)
    events = mne.read_events(event_fname)
    for evcnt in range(events.shape[0]-1):
        if events[evcnt,2] in np.array([2]):
            if events[evcnt+1,2] in np.array([79,90]):
                events[evcnt,2]=230
        if events[evcnt,2] in np.array([1]):
             if events[evcnt+1,2] in np.array([71]):
                 events[evcnt,2]=230
#    stim_delay = 0.034 # delay in s
#    events[:,0] = events[:,0]+np.round( raw.info['sfreq']*stim_delay )
    sfreq=raw.info['sfreq']/1000.
    sonset_raw=int(abs(tmin*1000)/sfreq)
    ###############################################################################
    # Read epochs for all channels, removing a bad one
    picks = mne.pick_types(raw.info, eeg=True, meg=True, eog=False, exclude='bads')
    print "Epoching"
    if subject_no in np.array([2,3,5,18]):#np.array([4,8,9]):
        reject = dict(eeg=150e-6, grad=200e-12, mag=4e-12)
    else:
        reject = dict(eeg=120e-6, grad=200e-12, mag=4e-12)#, eog=150e-6)
    event_ids = {'visual': 1, 'hear': 2, 'hand': 3}#,'pwordc': 6}#, 'neutral': 4, 'emotional': 5,'pwordc': 6}
    epochs = mne.Epochs(raw, events, event_ids, tmin, tmax, picks=picks, proj=True, baseline=(None, 0), reject=reject)
#    elenorig=len(epochs)
    epochs.drop_bad_epochs()
    elengood=len(epochs)
    print elengood/(450.-9.)
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
    epochs_list=range(nevents)
    for eid,event_id in zip(range(nevents),np.array([1,2,3])):
#        event_id = eid+1  
        epochs_list[eid]=copy.deepcopy(epochs[epochs.events[:,2]==event_id])
    equalize_epoch_counts(epochs_list,method='mintime')
    print 3*len(epochs_list[2])/(450.-9.)
    conind=0
    for event_no in which_conds:
#        conind=conind+1
#        print conind
        this_evoked = epochs_list[event_no].average()
#        epdata=this_evoked.data 
#        em=np.abs(epdata)
#        #emax=np.mean(em[:,:,600:900],axis=2)#
#        sa=em[:,int(sonset_raw+50/sfreq):int(sonset_raw+450/sfreq)]#.max(axis=2)
#        #bm=np.sqrt(np.mean(epdata[:,:,100:500]**2,axis=2))
#        bv=np.var(em[:,int(sonset_raw-250/sfreq):sonset_raw],axis=1)
#        edv=sa.copy()
#        for ii in range(sa.shape[1]):
#            edv[:,ii]=(sa[:,ii]**2)/bv    
#        ediv=np.sqrt(edv)
        snr = 3.0#np.mean(ediv)
        print snr
        lambda2 = 1.0 / snr ** 2
        method="MNE"
        fname_fwd = data_path + meg + 'ico5_forward_5-3L-EMEG-fwd.fif'
        forward = mne.read_forward_solution(fname_fwd,surf_ori=True, force_fixed=True)
        print "load source estimates"
        method = "MNE"  # use dSPM method (could also be MNE or sLORETA)
        fname_inv = data_path + meg + 'InvOp_ico5newreg_filtnew_pnt1_30_clean_ica_EMEG-inv.fif'
        # Load data
        inverse_operator = read_inverse_operator(fname_inv)
        srcin='/imaging/rf02/TypLexMEG/fsaverage/bem/fsaverage-ico-5-src.fif'
        src = mne.read_source_spaces(srcin)#inverse_operator['src']  # the source space used
#        fname = data_path + meg + 'firstMorphed_ico_SemDec_pnt1_30ica_'+event_names[event_no]+'_Source_Evoked_m300_600_5ms'  
        this_condition_normal = apply_inverse(this_evoked, inverse_operator, lambda2, method,pick_ori="normal")#mne.read_source_estimate(fname)
#        this_condition_none = apply_inverse(this_evoked, inverse_operator, lambda2, method,pick_ori=None)#mne.read_source_estimate(fname)

#        ntimes=this_condition.data.shape[1]
        print "extract ROI timecourse"
        subject_from=subjects[subject_no]
        subject_to='fsaverage'
        vertices_to = [np.arange(10242), np.arange(10242)]
        this_morphmat=mne.compute_morph_matrix(subject_from, subject_to, this_condition_normal.vertices, vertices_to, subjects_dir=data_path)
        this_stc_normal=this_condition_normal.morph_precomputed(subject_to,vertices_to,this_morphmat, subject_from)
#        this_stc_none=this_condition_none.morph_precomputed(subject_to,vertices_to,this_morphmat, subject_from)
        
        label_ts_normal = mne.extract_label_time_course(this_stc_normal, labellist, src, mode='mean_flip',return_generator=False)
#        label_ts_none = mne.extract_label_time_course(this_stc_none, labellist, src, mode='mean',return_generator=False)

        tcmat_normf[index,:,:,event_no]=label_ts_normal.copy()
#        tcmat_nonm[index,:,:,event_no]=label_ts_none.copy()
        
this_pathto=out_path+'wbw_newfilt_ROItcs_normal_meanflip_to30_ROI_ATAGIFGMTGVISHRHN_evVisHrHn.mat'#'DKA68_orig_modified_coh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot            
scipy.io.savemat(this_pathto,{'tcmat_normf':tcmat_normf})
#ltc=scipy.io.loadmat(this_pathto,mat_dtype=True); tcmat_normf=np.squeeze(ltc['tcmat_normf'])            

#this_pathto=out_path+'ROItcs_none_mean_ROI_ATAGIFGMTGVISHRHN_evVisHrHn.mat'#'DKA68_orig_modified_coh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot            
#scipy.io.savemat(this_pathto,{'tcmat_nonm':tcmat_nonm})
tcmat_normf_all=np.zeros((n_subjects,ntimes,nrois*nconds))
for index, meg in enumerate(ll):
    for event_no in which_conds:
        tcmat_normf_all[index,:,nrois*event_no:nrois*(event_no+1)]=-tcmat_normf[index,:,:,event_no].copy().T

plt.figure()

roinos=[4,5,6]
connos1=[0,1,2]
connos2=[0,1,2]

#roino=4
#conno1=0
#conno2=1
tlowlim=np.arange(350,751,25)
thilim=np.arange(375,776,25)
for roino in roinos:
    print roino
    for conno1 in connos1:
        print conno1
        for conno2 in np.setdiff1d(np.asarray(connos2),np.asarray([conno1])):
            print conno2
            X=np.zeros((n_subjects,len(tlowlim)))
            #            tst=350
            Xall=np.zeros((n_subjects,len(tlowlim),21))
            ii=-1
            for tl,th in zip (tlowlim,thilim):
                ii=ii+1
                X[:,ii]=np.mean(tcmat_normf[:,roino,tl:th,conno1],1)-np.mean(tcmat_normf[:,roino,tl:th,conno2],1)
                Xall[:,ii,:]=np.mean(tcmat_normf_all[:,tl:th,:],1)
                #                tst=tst+25
            #plt.plot(np.mean(X,0))
            #Xall=np.zeros(19,16,21)
            #            t=mne.stats.ttest_1samp_no_p(X)#
            #            t1,p,h=mne.stats.permutation_t_test(X)
            #plt.plot(p)
            #threshold=2
            #t, c, cp, h0=mne.stats.permutation_cluster_1samp_test(X, threshold=None, n_permutations=5000, tail=0, connectivity=None, verbose=None, n_jobs=1, seed=None, max_step=1, exclude=None, step_down_p=0, t_power=1, out_type='mask', check_disjoint=False, buffer_size=1000)    
            #
            #connectivity = spatial_tris_connectivity(grade_to_tris(5))   
            #connectivity = spatio_temporal_tris_connectivity(grade_to_tris(5), n_times=5)   
            n_permutations=10000
            
            cstep=0.1#(np.abs(X).max()-np.abs(X).min())/100
            cstart=0#np.abs(X).min()
            #p_threshold = 0.043
#            from scipy import stats as stats
#            p_threshold=0.05
#            t_threshold = -stats.distributions.t.ppf(p_threshold/2., n_subjects - 1)#dict(start=0, step=.1)#
#            #t_threshold=2
#            #            print('Clustering.')
#            tail=0
#            max_step=5;
#            T_obs, clusters, cluster_p_values, H0 = mne.stats.permutation_cluster_1samp_test(X, connectivity=None,  threshold=t_threshold, tail=tail, t_power=1, step_down_p=0.05, n_permutations=n_permutations, n_jobs=4)#, max_step=0)#, max_step=0)#spatial_exclude=exclude_ver, 
#            print cluster_p_values
#            print clusters
#plt.plot(T_obs)

#        method = 'MNE'  # can be 'MNE' or 'sLORETA'
#        mode = 'svd'
#        n_svd_comp = 1
#        labels = [labellist[jj].morph(subject_from='fsaverage', subject_to=subjects[subject_no], subjects_dir=data_path) for jj in range(len(labellist))]
#        stc_ctf_mne=psf_ctf_28Oct15.cross_talk_function(inverse_operator, forward, labels,method=method, lambda2=lambda2, signed=False, mode=mode, n_svd_comp=1, verbose=None)
#        
#        subject_from = subjects[subject_no]
#        subject_to = 'fsaverage'
#        vertices_to = [np.arange(10242), np.arange(10242)]
#        morphmat1=mne.compute_morph_matrix(subject_from, subject_to, stc_ctf_mne.vertices, vertices_to, subjects_dir=data_path)
#        stc_to=stc_ctf_mne.morph_precomputed(subject_to,vertices_to,morphmat1, subject_from)
#        Matx=stc_to.data[:,:-1]# Nv x Nlabels+1
#        labels_tc=np.zeros((len(labellist),n_times))
#        for this_labelc, this_label in enumerate(labellist):            
#            winner_verts=np.where(Matx[:,this_labelc]>=np.max(Matx[:,this_labelc])/2.)[0]
#            thisscale=1#Matx[winner_verts,this_labelc]/np.max(Matx[winner_verts,this_labelc])
#            label_predata=this_stc.data[winner_verts,:]
#            label_data1=np.tile(thisscale,[ntimes,1]).transpose()*label_predata
#            active_ones=np.where(np.mean(label_data1[:,sonset:],1)>=np.mean(np.mean(label_data1[:,sonset:],1)))[0]
#            label_data=label_data1[active_ones,:]
#            ini_tc=np.mean(label_data[:,:], axis=0)
#            labels_tc[this_labelc]=(ini_tc-np.min(ini_tc))/(np.max(ini_tc)-np.min(ini_tc))#ini_tc#zscore(ini_tc)#ini_tc#(ini_tc-np.min(ini_tc))/(np.max(ini_tc)-np.min(ini_tc))
#            X[index,:,conind]=(ini_tc-np.min(ini_tc))/(np.max(ini_tc)-np.min(ini_tc))#ini_tc#zscore(ini_tc)#ini_tc#(ini_tc-np.min(ini_tc))/(np.max(ini_tc)-np.min(ini_tc))
#            conind=conind+1
#cl_ideal[cntcl]=simstc_morphed_ideal[cntlts].in_label(labels[cntcl]).data
#clw_ideal[cntcl]=np.argmax(np.mean(np.abs(cl_ideal[cntcl][:,125:]),axis=1))
#thislabelsign2=mne.label.label_sign_flip(labels[cntcl],src_avg)[clw_ideal[cntcl]]
#label_ts_ideal[cntlts][cntcl,:]=thislabelsign2*np.expand_dims(cl_ideal[cntcl][clw_ideal[cntcl],:].transpose(),0)
#ATL 0, 7, 14
#AG 1, 8, 15
#IFG: 2, 9, 16
#MTG: 3, 10, 17
#vis: 4, 11, 18
#aud: 5, 12, 19
#hnd: 6, 13, 20
X_list=[np.squeeze(x) for x in np.split(Xall[:,:,np.array([4,6,18,20])], 4, axis=-1)]#np.array([1,3])#ROI3 in cond 1, 2, 3
#0-6,7-13,14-20
factor_levels = [2,2]  # number of levels in each factor
effects = 'A:B'  # this is the default signature for computing all effects
return_pvals = True  
pthresh = 0.05
n_permutations=10000
##    max_step=1;
f_thresh = mne.stats.f_threshold_mway_rm(n_subjects, factor_levels, effects, pthresh)
thresh_tfce=dict(start=0, step=0.05)
#        
def stat_fun(*args):  
    return f_mway_rm(np.swapaxes(args, 1, 0), factor_levels=factor_levels,effects=effects, return_pvals=return_pvals)[0]
T_obs, clusters, cluster_p_values, H0 = mne.stats.permutation_cluster_test(X_list, connectivity=None,  threshold=f_thresh, stat_fun=stat_fun,  t_power=1,  n_permutations=n_permutations, n_jobs=4,max_step=0)#,step_down_p=0.05, max_step=0)#, max_step=0)#spatial_exclude=exclude_ver, 
plt.figure()
plt.plot(tlowlim-300,T_obs,color='black')
sigline=np.zeros((len(tlowlim),))
sigline[:]=np.NAN
sigline[7:8]=0.5
plt.plot(tlowlim-300,sigline,'red')
plt.xlabel('time (ms)')
plt.ylabel('F value')
print "interaction"
print clusters
print cluster_p_values

thisX=np.swapaxes(Xall[:,:,np.array([4,6,18,20])],1,2)#subxcon
f, p= f_mway_rm(thisX, factor_levels=factor_levels,effects=effects, return_pvals=return_pvals)
h, pc=mne.stats.fdr_correction(p)

#source_space = grade_to_tris(5)
## as we only have one hemisphere we need only need half the connectivity
#print('Computing connectivity.')
#connectivity = spatial_tris_connectivity(source_space)
#
##    Now let's actually do the clustering. Please relax, on a small
##    notebook and one single thread only this will take a couple of minutes ...
#pthresh = 0.05
###    max_step=1;
#f_thresh = mne.stats.f_threshold_mway_rm(n_subjects, factor_levels, effects, pthresh)
##
##    To speed things up a bit we will ...
#n_permutations = 5000
#fsave_vertices = [np.arange(10242), np.arange(10242)]
#T_obs, clusters, cluster_p_values, H0 =  permutation_cluster_test(X_list, threshold=f_thresh, tail=1,stat_fun=stat_fun, n_permutations=n_permutations, n_jobs=4)
#p_thr=0.05
#good_cluster_inds = np.where(cluster_p_values <p_thr)[0]
#print cluster_p_values[good_cluster_inds]; print good_cluster_inds
#
#
#
#
#xx1=np.mean(X,0)
#cond1=xx1[:,:5]
#cond2=xx1[:,5:10]
#roi1=xx1[:,np.array([4,9,14])]
#roi2=xx1[:,np.array([2,7,12])]
#plt.plot(roi1)
#plt.figure(2)
#plt.plot(roi2)

#    Matx_zscore=np.zeros(Matx.shape)
##    for cntr in range(Matx.shape[1]-1):
##        Matx_normalised[:,cntr]=np.divide(Matx[:,cntr],Matx[:,-1])
#    Matx_zscore=zscore(Matx,axis=1)
#    assign_vertex=np.argmax(Matx_zscore, axis=1)
#    msort=np.sort(Matx_zscore,axis=1)
#    mtokeep=msort[:,-1]-msort[:,-2]
#    for avii in range(len(assign_vertex)):
#        if Matx[avii,assign_vertex[avii]]<np.max(Matx[:,assign_vertex[avii]])/4.:
#            assign_vertex[avii]=-100
#            
#    assign_vertex[np.logical_or(msort[:,-1]<3, mtokeep<=1)]=-100
#    print "ctf finished!"
#    
#    
#    
## pick MEG channels
#    picks = mne.pick_types(raw.info, meg=True, eeg=True, stim=False, eog=True, exclude='bads')
#
## Read epochs
#    epochs = mne.Epochs(raw, events, event_ids, tmin, tmax, picks=picks, baseline=(None, 0), proj=True, reject=reject)
#    epochs_concrete = epochs['cncrt_wrd']
#    epochs_abstract = epochs['abs_wrd']
#
#    mne.epochs.equalize_epoch_counts([epochs_concrete,epochs_abstract])
#    epochs=epochs_abstract
#    events_no=epochs.copy().get_data().shape[0]
#    stcs=range(events_no)
#    for ccc in range(len(stcs)):
#        print ccc		
#        data_in = data_path + meg + 'Morphed_ico_SemLoc_epochs' + str(ccc) + '_icaclean_Abstract_Source_Evoked_m500_700'
#        stc_orig=mne.read_source_estimate(data_in)
#        #stc_orig.resample(100)
#        stcs[ccc]=stc_orig
#
## First, we find the most active vertex in the left auditory cortex, which
##  we will later use as seed for the connectivity computation
##    snr = 3.0
##    lambda2 = 1.0 / snr ** 2
#    evoked = epochs.average()
#    stc_incoh = data_path + meg + 'firstMorphed_ico_SemLoc_ica_Abstract_Source_Evoked_m500_700'
#    stc_coh = mne.read_source_estimate(stc_incoh)
#    #stc_coh.resample(100)
#
#    for label1,lii in zip(labellist,label_index):
#        print label1
#        tmin, tmax = -0.5, 0.7
#        if subject_inds[0] in np.array([3,4,13]):
#            reject = dict(eeg=150e-6, grad=200e-12, mag=4e-12) #eog=150e-6,
#        else:
#            reject = dict(eeg=120e-6, grad=200e-12, mag=4e-12) #eog=150e-6,
#        event_ids = {'cncrt_wrd': 1, 'abs_wrd': 2}
## frequency bands with variable number of cycles for wavelets
#        frequencies = np.arange(8, 45, 1)
#        n_cycles = frequencies / float(7)
#        n_cycles[frequencies<=20] = 2
#        n_cycles[frequencies<=20] += np.arange(0,1,0.077)
## n_cycles[np.logical_and((frequencies>10), (frequencies<=20))] = 3
#        subject_from = subjects[subject_inds[0]]
#        subject_to = 'fsaverage'
#        vertices_to = [np.arange(10242), np.arange(10242)]
#
#        method = "MNE"  # use dSPM method (could also be MNE or sLORETA)
#
#        #fname_label = label_path +  label_name + '.label'
#        #label1 = mne.read_label(fname_label)
#        actual_label=stc_to.in_label(label1)
#        if label1.name.endswith('-lh'):
#            actual_label_vert=actual_label.vertices[0]
#        else:
#            actual_label_vert=actual_label.vertices[1]+len(stc_to.lh_vertno)
#            
#        best_verts=np.intersect1d(np.where(assign_vertex==lii)[0],actual_label_vert)
#        bestvert_values=Matx_zscore[best_verts,lii]
#        best_101thvert=sorted(bestvert_values,reverse=True)[1]
#        best_100vert=best_verts#best_verts[np.where(bestvert_values>best_101thvert)]
#        print Matx_zscore[best_100vert,lii].mean()
#        label_stc_data=np.zeros((stc_to.data.shape[0],1))
#        label_stc_data[best_100vert]=1.
#        label_stc=mne.SourceEstimate(label_stc_data, vertices=vertices_to,tmin=1e-3 * 1., tstep=1e-3 * 1., subject='fsaverage')
#        label_new=mne.stc_to_label(label_stc, smooth=False,subjects_dir=data_path,connected=False)[0]
#        label_new.name=label1.name[:-3]+'_'+str(subject_no)+'_ctf-lh'
#        label_new_path=label_path+label_new.name
#        #label_new.save(label_new_path)
#        print label_new.vertices.shape; print "label tc extracted"
#        seed_ts = mne.extract_label_time_course(stcs, label_new, src, mode='pca_flip')
#        comb_ts = zip(seed_ts, stcs)
#        
## Restrict the source estimate to the label in the left auditory cortex
#        #stcp = apply_inverse(evoked, inverse_operator, lambda2, method, pick_ori="normal")
#        #stc = mne.morph_data(subject_from, subject_to, stcp, subjects_dir=data_path, n_jobs=1, grade=vertices_to)
##        stc_label = stc_coh.in_label(label1)
##
### Find number and index of vertex with most power
##        src_pow = np.sum(stc_label.data ** 2, axis=1)
##        if label_name[0] is 'l':
##            print label_name
##            seed_vertno = stc_label.vertices[0][np.argmax(src_pow)]
##            seed_idx = np.searchsorted(stc_coh.vertices[0], seed_vertno)  # index in original stc
##            print seed_idx
##
##        else:
##            print label_name
##            seed_vertno = stc_label.vertices[1][np.argmax(src_pow)]
##            seed_idx = stc_coh.vertices[0].shape[0]+ np.searchsorted(stc_coh.vertices[1], seed_vertno)  # index in original stc
##            print seed_idx
#
## Generate index parameter for seed-based connectivity analysis
##        n_sources = stc_coh.data.shape[0]
##        indices = seed_target_indices([0], np.arange(n_sources))
#        vertices = [src[i]['vertno'] for i in range(2)]
#        n_signals_tot = 1 + len(vertices[0]) + len(vertices[1])
#
#        indices = seed_target_indices([0], np.arange(1,n_signals_tot))
#
#
## Compute inverse solution and for each epoch. By using "return_generator=True"
## stcs will be a generator object instead of a list. This allows us so to
## compute the Coherence without having to keep all source estimates in memory.
#
## Now we are ready to compute the Coherence in the alpha and beta band.
## fmin and fmax specify the lower and upper freq. for each band, resp.
#        fmin = (4., 8., 13., 31.)
#        fmax = (7., 12., 30, 45.)
#        tmin2 = 0.0
#        tmax2 = 0.551
#        tminb2 = -0.4
#        tmaxb2 = 0.
#        sfreq = raw.info['sfreq']  # the sampling frequency
#        #cwt_n_cycles=np.array([2,2,2,2,3,3,3,3,3,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7])
#        cwt_frequencies = np.arange(4, 46, 1)#np.array([6,10,22,38])#
#        cwt_n_cycles = cwt_frequencies / float(3)#np.array([2,3,5,7])#
#        #cwt_n_cycles[cwt_frequencies<=20] = 2
#        #cwt_n_cycles[cwt_frequencies<=20] += np.arange(0,1,0.059)#0.117)#0.059)
#        #con_blcn, freqs_blcn, times_blcn, n_epochs_blcn, n_tapers_blcn = spectral_connectivity(label_ts_concrete, method=con_methods, mode='multitaper', sfreq=sfreq, fmin=fmin, fmax=fmax, faverage=True, tmin=tminb2, tmax=tmaxb2, mt_adaptive=True, n_jobs=2)
#    
#        #con_blab, freqs_bl, times_bl, n_epochs_bl, n_tapers_bl = spectral_connectivity(label_ts_abstract, method=con_methods, mode='multitaper', sfreq=sfreq, fmin=fmin, fmax=fmax, faverage=True, tmin=tminb2, tmax=tmaxb2, mt_adaptive=True, n_jobs=2)
#    
#
## Now we compute connectivity. To speed things up, we use 2 parallel jobs
## and use mode='multitaper', mt_adaptive=True, which uses a FFT with a Hanning window
## to compute the spectra (instead of multitaper estimation, which has a
## lower variance but is slower). By using faverage=True, we directly
## average the Coherence in the alpha and beta band, i.e., we will only
## get 2 frequency bins
#	
#        print "spectral connectivity coherence"
#        coh=np.zeros((indices[0].shape[0],4,2))
#        for actc,tc in enumerate(range(150,251,100)): #time count
#            tminw=tc/1000.
#            tmaxw=(tc+200)/1000.
#            coh[:,:,actc], freqs1, times1, n_epochs1, _ = spectral_connectivity(comb_ts, method='coh', mode='multitaper', indices=indices,sfreq=sfreq, fmin=fmin, fmax=fmax, faverage=True, tmin=tminw, tmax=tmaxw, n_jobs=4)
###spectral_connectivity(data, method='coh', indices=None, sfreq=6.283185307179586, mode='multitaper', fmin=None, fmax=inf, fskip=0, faverage=False, tmin=None, tmax=None, mt_bandwidth=None, mt_adaptive=False, mt_low_bias=True, cwt_frequencies=None, cwt_n_cycles=7, block_size=1000, n_jobs=1, verbose=None)
#        print('coh fin ')
#
## Generate a SourceEstimate with the Coherence. This is simple since we
## used a single seed. For more than one seeds we would have to split coh.
## Note: We use a hack to save the frequency axis as time
#        print "use a hack to save the frequency axis as time"
#        tmin1 = 0.150#np.mean(freqs[0])
#        tstep1 = 0.100#np.mean(freqs[1]) - tmin
#        coh_stc = mne.SourceEstimate(coh[:,0,:], vertices=stc_coh.vertices, tmin= tmin1, tstep= tstep1, subject='fsaverage')
#        data_out = data_path + meg +  'Morphed_ico_SemLoc_ica_equalized_Abstract_Mlttap_PCA_ctf_Coherence_Theta_150_450_200ms_' + label1.name[:-3]
#        print 'theta'
#        coh_stc.save(data_out)
#        
#        coh_stc = mne.SourceEstimate(coh[:,1,:], vertices=stc_coh.vertices, tmin= tmin1, tstep= tstep1, subject='fsaverage')
#        data_out = data_path + meg +  'Morphed_ico_SemLoc_ica_equalized_Abstract_Mlttap_PCA_ctf_Coherence_Alpha_150_450_200ms_' + label1.name[:-3]
#        print 'alpha'
#        coh_stc.save(data_out)
#        
#        coh_stc = mne.SourceEstimate(coh[:,2,:], vertices=stc_coh.vertices, tmin= tmin1, tstep= tstep1, subject='fsaverage')
#        data_out = data_path + meg +  'Morphed_ico_SemLoc_ica_equalized_Abstract_Mlttap_PCA_ctf_Coherence_Beta_150_450_200ms_' + label1.name[:-3]
#        print 'beta'
#        coh_stc.save(data_out)
#        
#        coh_stc = mne.SourceEstimate(coh[:,3,:], vertices=stc_coh.vertices, tmin= tmin1, tstep= tstep1, subject='fsaverage')
#        data_out = data_path + meg +  'Morphed_ico_SemLoc_ica_equalized_Abstract_Mlttap_PCA_ctf_Coherence_Gamma_150_450_200ms_' + label1.name[:-3]
#        print 'gamma'
#        coh_stc.save(data_out)
#	
##############
#
#        print "spectral connectivity ppc"
#        dw_pli=np.zeros((indices[0].shape[0],4,2))
#        for actc,tc in enumerate(range(150,251,100)): #time count
#            tminw=tc/1000.
#            tmaxw=(tc+200)/1000.
#            dw_pli[:,:,actc], freqs1, times1, n_epochs1, _ = spectral_connectivity(comb_ts, method='wpli2_debiased', mode='multitaper', indices=indices,sfreq=sfreq, fmin=fmin, fmax=fmax, faverage=True, tmin=tminw, tmax=tmaxw, n_jobs=4)
###spectral_connectivity(data, method='coh', indices=None, sfreq=6.283185307179586, mode='multitaper', fmin=None, fmax=inf, fskip=0, faverage=False, tmin=None, tmax=None, mt_bandwidth=None, mt_adaptive=False, mt_low_bias=True, cwt_frequencies=None, cwt_n_cycles=7, block_size=1000, n_jobs=1, verbose=None)
#        print('ppc fin ')
#
## Generate a SourceEstimate with the Coherence. This is simple since we
## used a single seed. For more than one seeds we would have to split coh.
## Note: We use a hack to save the frequency axis as time
#        print "use a hack to save the frequency axis as time"
#        tmin1 = 0.150#np.mean(freqs[0])
#        tstep1 = 0.100#np.mean(freqs[1]) - tmin
#        coh_stc = mne.SourceEstimate(dw_pli[:,0,:], vertices=stc_coh.vertices, tmin= tmin1, tstep= tstep1, subject='fsaverage')
#        data_out = data_path + meg +  'Morphed_ico_SemLoc_ica_equalized_Abstract_Mlttap_PCA_ctf_WPLI2d_Theta_150_450_200ms_' + label1.name[:-3]
#        print 'theta'
#        coh_stc.save(data_out)
#        
#        coh_stc = mne.SourceEstimate(dw_pli[:,1,:], vertices=stc_coh.vertices, tmin= tmin1, tstep= tstep1, subject='fsaverage')
#        data_out = data_path + meg +  'Morphed_ico_SemLoc_ica_equalized_Abstract_Mlttap_PCA_ctf_WPLI2d_Alpha_150_450_200ms_' + label1.name[:-3]
#        print 'alpha'
#        coh_stc.save(data_out)
#        
#        coh_stc = mne.SourceEstimate(dw_pli[:,2,:], vertices=stc_coh.vertices, tmin= tmin1, tstep= tstep1, subject='fsaverage')
#        data_out = data_path + meg +  'Morphed_ico_SemLoc_ica_equalized_Abstract_Mlttap_PCA_ctf_WPLI2d_Beta_150_450_200ms_' + label1.name[:-3]
#        print 'beta'
#        coh_stc.save(data_out)
#        
#        coh_stc = mne.SourceEstimate(dw_pli[:,3,:], vertices=stc_coh.vertices, tmin= tmin1, tstep= tstep1, subject='fsaverage')
#        data_out = data_path + meg +  'Morphed_ico_SemLoc_ica_equalized_Abstract_Mlttap_PCA_ctf_WPLI2d_Gamma_150_450_200ms_' + label1.name[:-3]
#        print 'gamma'
#        coh_stc.save(data_out)

###############
#        print "spectral connectivity psi"
#        psi=np.zeros((n_sources,4,2))
#        for actc,tc in enumerate(range(150,251,100)): #time count
#            tminw=tc/1000.
#            tmaxw=(tc+200)/1000.
#            psi[:,:,actc], freqs1, times1, n_epochs1, n_tapers1 = mne.connectivity.phase_slope_index(comb_ts, mode='multitaper', mt_adaptive=True, indices=indices,sfreq=sfreq, fmin=fmin, fmax=fmax,  tmin=tminw, tmax=tmaxw, n_jobs=4)
###spectral_connectivity(data, method='coh', indices=None, sfreq=6.283185307179586, mode='multitaper', fmin=None, fmax=inf, fskip=0, faverage=False, tmin=None, tmax=None, mt_bandwidth=None, mt_adaptive=False, mt_low_bias=True, cwt_frequencies=None, cwt_n_cycles=7, block_size=1000, n_jobs=1, verbose=None)
#        print('dwpli fin ')
#
## Generate a SourceEstimate with the Coherence. This is simple since we
## used a single seed. For more than one seeds we would have to split coh.
## Note: We use a hack to save the frequency axis as time
#        print "use a hack to save the frequency axis as time"
#        tmin1 = 0.150#np.mean(freqs[0])
#        tstep1 = 0.100#np.mean(freqs[1]) - tmin
#        coh_stc = mne.SourceEstimate(psi[:,0,:], vertices=stc_coh.vertices, tmin= tmin1, tstep= tstep1, subject='fsaverage')
#        data_out = data_path + meg +  'firstMorphed_SemLoc_icaclean_equalized_Concrete_Multitaper_PSI_Theta_150_450_200ms_' + label_name[0:-3]
#        print 'theta'
#        coh_stc.save(data_out)
#        
#        coh_stc = mne.SourceEstimate(psi[:,1,:], vertices=stc_coh.vertices, tmin= tmin1, tstep= tstep1, subject='fsaverage')
#        data_out = data_path + meg +  'firstMorphed_SemLoc_icaclean_equalized_Concrete_Multitaper_PSI_Alpha_150_450_200ms_' + label_name[0:-3]
#        print 'alpha'
#        coh_stc.save(data_out)
#        
#        coh_stc = mne.SourceEstimate(psi[:,2,:], vertices=stc_coh.vertices, tmin= tmin1, tstep= tstep1, subject='fsaverage')
#        data_out = data_path + meg +  'firstMorphed_SemLoc_icaclean_equalized_Concrete_Multitaper_PSI_Beta_150_450_200ms_' + label_name[0:-3]
#        print 'beta'
#        coh_stc.save(data_out)
#        
#        coh_stc = mne.SourceEstimate(psi[:,3,:], vertices=stc_coh.vertices, tmin= tmin1, tstep= tstep1, subject='fsaverage')
#        data_out = data_path + meg +  'firstMorphed_SemLoc_icaclean_equalized_Concrete_Multitaper_PSI_Gamma_150_450_200ms_' + label_name[0:-3]
#        print 'gamma'
#        coh_stc.save(data_out)
#
