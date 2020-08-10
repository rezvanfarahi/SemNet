# -*- coding: utf-8 -*-
"""
Created on Wed Dec 23 14:34:39 2015

@author: rf02
"""

"""
==========================================================
Compute point-spread functions (PSFs) for MNE/dSPM/sLORETA
==========================================================

PSFs are computed for four labels in the MNE sample data set
for linear inverse operators (MNE, dSPM, sLORETA).
PSFs describe the spread of activation from one label
across the cortical surface.
"""

# Authors: Olaf Hauk <olaf.hauk@mrc-cbu.cam.ac.uk>
#          Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
#
# License: BSD (3-clause)



# Author: Martin Luessi <mluessi@nmr.mgh.harvard.edu>
#
# License: BSD (3-clause)

print(__doc__)
print "hi! this script is for spectral functional connectivity. Very nice maps are fruits which worth waiting for I bet!"
import sys

#sys.path.append('/imaging/local/software/python_packages/nibabel/1.3.0')
#sys.path.append('/imaging/local/software/python_packages/pysurfer/v0.4')
#sys.path.append('/imaging/rf02/scikit-learn-0.15.0')
#
##sys.path.append('/imaging/local/software/python_packages/scikit-learn/v0.14.1')
#sys.path.append('/imaging/local/software/mne_python/la_v0.9')
sys.path.insert(1,'/imaging/rf02/mne_python_11')
#sys.path.insert(1,'/imaging/local/software/mne_python/la_v0.9full')
#sys.path.insert(1,'/imaging/local/software/mne_python/v0.11')


sys.path.insert(1,'/imaging/rf02/scikit-learn-0.16.1')
sys.path.insert(1,'/home/rf02/.local/lib/python2.7/site-packages')


# for qsub
# (add ! here if needed) 
sys.path.append('/imaging/local/software/EPD/latest/x86_64/bin/python')
#sys.path.append('/imaging/rf02/TypLexMEG/mne_python_v8/')
#sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
###
import sklearn
from scipy.stats.mstats import zscore
import os
import numpy as np
import mne
from mne import io
from mne.io import Raw
import operator
from mne.minimum_norm import (read_inverse_operator, make_inverse_operator, write_inverse_operator, apply_inverse, apply_inverse_epochs)
from mne.minimum_norm.inverse import (prepare_inverse_operator)
import scipy.io
from mne import find_events
from mne.connectivity import seed_target_indices, spectral_connectivity
from mne import read_labels_from_annot
from mne.label import _n_colors as nc
from mne.label import _read_annot as mnera
from mne.label import _write_annot as mnewa
from mne.viz import circular_layout, plot_connectivity_circle
import matplotlib.pyplot as plt
import copy
import pickle
print mne.__path__
def remove_fake_R(labels, Rmat, con_mat, psd_mat ):
    for cnt1 in range(len(labels)): #seed label
        for cnt2 in range(len(labels)): #destination labels
            if cnt1<cnt2:
                con_mat[cnt1,cnt2,:]=con_mat[cnt2,cnt1,:]
    con_mat_mean=np.mean(con_mat,axis=2)
    sRmat=np.divide(Rmat,np.max(Rmat,axis=0))  
    label_patches=np.where(sRmat>0.25)   
    con_mat_modified=np.zeros((len(labels),len(labels),con_mat.shape[2]))        
    for cnt1 in range(len(labels)): #seed label
        print cnt1
        for cnt2 in range(len(labels)): #destination labels
            if cnt1!=cnt2:
                this_patch=label_patches[0][label_patches[1]==cnt2]
                thislabw=np.where(this_patch==cnt2)[0]
                notthislabw=np.where(this_patch!=cnt2)[0]
                coef_mat_all=Rmat[this_patch,:][:,this_patch].transpose()
                coef_mat=coef_mat_all[notthislabw,:][:,notthislabw]
                this_coef_vect=coef_mat_all[thislabw,notthislabw]
                est_corrs=np.zeros((psd_mat.shape[1],1))
                for fcnt in range(psd_mat.shape[1]):
                    Gzz=psd_mat[cnt1,fcnt]*np.ones((len(this_patch)-1,))
                    scale_coeff0=np.sqrt(Gzz*psd_mat[this_patch[notthislabw],fcnt])
                    psd_patch=psd_mat[this_patch[notthislabw],fcnt]
                    scale_coeff=np.divide(coef_mat,scale_coeff0[:,np.newaxis])
                    con_vect=con_mat[this_patch[notthislabw],cnt1,fcnt]
                    est_Sxys=np.divide(con_vect,np.diag(scale_coeff))#np.linalg.solve(scale_coeff,con_vect)
                    est_Sxy=np.max(est_Sxys*this_coef_vect)#np.sum(est_Sxys*this_coef_vect)
                    est_corrs[fcnt]=est_Sxy/np.sqrt(Gzz[0]*psd_mat[this_patch[thislabw],fcnt])
                    est_corr=est_corrs[fcnt]#np.mean(est_corrs)
                    act_corr=con_mat[cnt1,cnt2,fcnt]
                    if np.abs(act_corr-est_corr)/act_corr<0.2:
                        con_mat_modified[cnt1,cnt2,fcnt]=0#np.sqrt(psd_mat[cnt2,fcnt]*Gzz[0])#
                    else:
                        con_mat_modified[cnt1,cnt2,fcnt]=act_corr
    con_mat_final=np.sign(con_mat_modified)*np.minimum(np.abs(con_mat_modified),np.abs(np.transpose(con_mat_modified,(1,0,2))))
    #    for cnt1 in range(len(labels)): #seed label
    #        for cnt2 in range(len(labels)): #destination labels
    #            if cnt1<cnt2:
    #               con_mat_final[cnt1,cnt2]=0 
    return con_mat_final
def estimate_actual_R(labels, Rmat, con_mat, psd_mat ):
    for cnt1 in range(len(labels)): #seed label
        for cnt2 in range(len(labels)): #destination labels
            if cnt1<cnt2:
                con_mat[cnt1,cnt2,:]=con_mat[cnt2,cnt1,:]
    
    sRmat=np.divide(Rmat,np.max(Rmat,axis=0))  
    label_patches=np.where(sRmat>0.25)   
    con_mat_modified=np.zeros((len(labels),len(labels)))        
    for cnt1 in range(len(labels)): #seed label
        print cnt1
        for cnt2 in range(len(labels)): #destination labels
            if cnt1!=cnt2:
                this_patch=label_patches[0][label_patches[1]==cnt2]
                thislabw=np.where(this_patch==cnt2)[0]
                coef_mat=Rmat[this_patch,:][:,this_patch].transpose()
                est_corr=np.zeros((psd_mat.shape[1],1))
                for fcnt in range(psd_mat.shape[1]):
                    Gzz=psd_mat[cnt1,fcnt]*np.ones((len(this_patch),))
                    scale_coeff0=np.sqrt(Gzz*psd_mat[this_patch,fcnt])
                    psd_patch=psd_mat[this_patch,fcnt]
                    scale_coeff=np.divide(coef_mat,scale_coeff0[:,np.newaxis])
                    con_vect=con_mat[this_patch,cnt1,fcnt]
                    est_Sxys=np.linalg.solve(scale_coeff,con_vect)
                    est_psds=np.linalg.solve(coef_mat,psd_patch)
                    est_Sxy=est_Sxys[thislabw]
                    est_psd=est_psds[thislabw]
                    est_corr[fcnt]=est_Sxy/np.sqrt(psd_mat[cnt2,fcnt]*Gzz[0])
#                    if est_psd>0:
#                        est_corr[fcnt]=est_Sxy/np.sqrt(est_psd*Gzz[0])#np.sqrt(psd_mat[cnt2,fcnt]*Gzz[0])#
#                    else:
#                        est_corr[fcnt]=0
                con_mat_modified[cnt1,cnt2]=np.mean(est_corr)
    con_mat_final=(con_mat_modified+np.transpose(con_mat_modified))/2.
#    for cnt1 in range(len(labels)): #seed label
#        for cnt2 in range(len(labels)): #destination labels
#            if cnt1<cnt2:
#               con_mat_final[cnt1,cnt2]=0 
    return con_mat_final
                

data_path = '/imaging/rf02/TypLexMEG/' # root directory for your MEG data
os.chdir(data_path)
event_path = '/group/erp/data/olaf.hauk/MEG/TypLexM/data/MF22_new/'
#label_path = '/group/erp/data/olaf.hauk/MEG/TypLexM/labels/createdlabels_SL/'
inv_path = '/imaging/rf02/TypLexMEG/'
subjects_dir = data_path 
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
ii=-1

label_path='/imaging/rf02/TypLexMEG/parcellation/labels_and_CTFs/'
results_path='/imaging/rf02/parcellation/leakage_pattern/'
#labellist_path = [label_path+'rostfront_smallseed-lh.label',label_path+'midtemp_smallseed-lh.label',label_path+'latocci_smallseed-lh.label']#[ 'lh.posterior_fusiform','rh.posterior_fusiform', 'lh.lateraloccipital','rh.lateraloccipital','lh.precentral','rh.precentral', 'lh.ATL_latmed','rh.ATL_latmed']#labellist = [ 'leftATL-lh','rightATL-rh']#, 'medialtempright-rh','medialtempleft-lh']'lh.postcentral','rh.postcentral',
#labellist_path = [label_path+'IFG_smallseed-lh.label',label_path+'IPL_smallseed-lh.label']#[ 'lh.posterior_fusiform','rh.posterior_fusiform', 'lh.lateraloccipital','rh.lateraloccipital','lh.precentral','rh.precentral', 'lh.ATL_latmed','rh.ATL_latmed']#labellist = [ 'leftATL-lh','rightATL-rh']#, 'medialtempright-rh','medialtempleft-lh']'lh.postcentral','rh.postcentral',
labellist_path = [label_path+'LO_smallseed-lh.label',label_path+'RMF_smallseed-lh.label',label_path+'STS_smallseed-lh.label']#[ 'lh.posterior_fusiform','rh.posterior_fusiform', 'lh.lateraloccipital','rh.lateraloccipital','lh.precentral','rh.precentral', 'lh.ATL_latmed','rh.ATL_latmed']#labellist = [ 'leftATL-lh','rightATL-rh']#, 'medialtempright-rh','medialtempleft-lh']'lh.postcentral','rh.postcentral',
labellist=[mne.read_label(label) for label in labellist_path]
for lbl in labellist:
    lbl.values.fill(1.0)
#
which_parc=1 #names below
parc_name=['aparc',
'fsaverage_finalpostparc_smooth_stddiff_split_merge_halfmax_86_badremoved_74',
'aparc.a2009s',
'fsaverage_finalpostparc2009_smooth_stddiff_split_merge_halfmax_98_badremoved_74_cmatch',
'fsaverage_RG_finalpostparc_smooth_stddiff_split_merge_halfmax_82_badremoved_70_cmatch'
]##'fsaverage_finalpostparc2009_smooth_stddiff_split_merge_halfmax_98_badremoved_74_cmatch'#'fsaverage_RG_finalpostparc_smooth_stddiff_split_merge_halfmax_82_badremoved_70_cmatch'
#
labels_pre = mne.read_labels_from_annot('fsaverage', parc=parc_name[which_parc], subjects_dir=subjects_dir)
#labels_pre.remove(labels_pre[-1])
#    labels_pre.remove(labels_pre[-1])
Rmat_names=['ctf_ypos_aparc.mat',
'ctf_ypos_finalpostparc_smooth_stddiff_split_merge_halfmax_86_badremoved_74.mat',
'ctf_ypos_anat2009.mat',
'ctf_ypos_finalpostparc2009_smooth_stddiff_split_merge_halfmax_98_badremoved_74_cmatch.mat',
'ctf_ypos_RG_finalpostparc_smooth_stddiff_split_merge_halfmax_82_badremoved_70_cmatch.mat'
]
path_Rmat=results_path+Rmat_names[which_parc]
ltc=scipy.io.loadmat(path_Rmat,mat_dtype=True); ctf_to_ypos=np.squeeze(ltc['ctf_to_ypos'])

ypos_names=['ypos_aparc.mat',
'ypos_aparcmod_badremoved_74.mat',
'ypos_aparc2009.mat',
'ypos_aparc2009mod_badremoved_74.mat',
'ypos_RG_badremoved_70.mat'
]
path_ypos=results_path+ypos_names[which_parc]
ltc=scipy.io.loadmat(path_ypos,mat_dtype=True); ypos_sort=np.squeeze(ltc['ypos_sort'])
labels=[labels_pre[jj] for jj in ypos_sort]

for cnt, meg in enumerate(ll):
    ii=ii+1
    print cnt
    subjects_dir = data_path 
    fname_fwd = data_path + meg + 'ico5_forward_5-3L-EMEG-fwd.fif'
    fname_inv = inv_path + meg + 'InvOp_ico5newreg_fft_1_48_clean_ica_EMEG-inv.fif'
    vertices_avg = [np.arange(10242), np.arange(10242)]
    subject_no=subject_inds[cnt]
    raw_fname = data_path + meg + 'semloc_ssstf_fft_1_48_clean_ica_raw.fif' #clean_ssp_
    event_fname = event_path + meg + 'semloc_raw_ssstnew.txt'
    subjects_dir = data_path 
    
    
    raw = io.Raw(raw_fname)
    events1 = mne.read_events(event_fname)
    stim_delay = 0.034 # delay in s
    events1[:,0]= events1[:,0]+np.round( raw.info['sfreq']*stim_delay )
    
    # Read epochs for all channels, removing a bad one
    print "epochs"
    picks = mne.pick_types(raw.info, eeg=True, meg=True, eog=False, exclude=[])
    epochs1 = mne.Epochs(raw, events1, tmin=-0.5, tmax=0.7,event_id=1, picks=picks, baseline=(None, 0), proj=True, reject=None, preload=True)
    
    forward = mne.read_forward_solution(fname_fwd,surf_ori=True, force_fixed=True)#, surf_ori=True)
    # read inverse operators
    inverse_operator_eegmeg = read_inverse_operator(fname_inv)
    info = copy.deepcopy(inverse_operator_eegmeg['info'])
    info['sfreq'] = 1000.  # necessary
    info['projs'] = inverse_operator_eegmeg['projs']
    ## simulate activity in labels
    morphed_labels=[labellist[jj].morph(subject_from='fsaverage', smooth=5, subject_to=subjects[subject_no], subjects_dir=data_path) for jj in range(len(labellist))]

    src_sub=forward['src']
    leadfield=forward['sol']['data']
    
    
    n_epochs=50
    LO_phase=np.pi*(2*np.random.random(n_epochs,)-1)
    RMF_phase=np.pi*(2*np.random.random(n_epochs,)-1)
    stc_sim_ideal=range(n_epochs)
    for ecnt in range(n_epochs):
        t=np.linspace(-6*np.pi,6*np.pi,3020/5)
        
        LO=np.concatenate((np.zeros(125,),1e-9*np.sin(t+LO_phase[ecnt])),axis=0)#+LO_phase[ecnt]
        RMF=np.concatenate((np.zeros(125,),1e-9*np.sin(2*t+RMF_phase[ecnt])),axis=0)#+RMF_phase[ecnt]
        MTG=np.concatenate((np.zeros(125,),1e-9*np.sin(t+LO_phase[ecnt]+(np.pi/3.))),axis=0)+np.concatenate((np.zeros(125,),1e-9*np.sin(2*t+RMF_phase[ecnt]-(np.pi/3.))),axis=0)
        this_std=np.sqrt(np.mean(RMF[125:]**2)/9.)
        datamat=np.concatenate((LO[:,np.newaxis],RMF[:,np.newaxis],MTG[:,np.newaxis]),axis=1).transpose()
        tmin=-0.124
        tstep=0.001
        lvert=src_sub[0]['vertno']
        rvert=src_sub[1]['vertno']
        vertices_sub = [lvert, rvert]
        thislabel_data=np.zeros((len(lvert)+len(rvert),len(LO)))
    #    thislabel_stc = mne.SourceEstimate(thislabel_data, vertices=vertices_sub, tmin=tmin, tstep=tstep, subject=subjects[subject_no], verbose=None)
    #    LOdat=thislabel_stc.in_label(morphed_labels[0])
    #    RMFdat=thislabel_stc.in_label(morphed_labels[1])
    #    MTGdat=thislabel_stc.in_label(morphed_labels[2])
      
        for llcnt,ll in enumerate(morphed_labels):
            if ll.hemi == 'rh':
                # for RH labels, add number of LH vertices
                offset = forward['src'][0]['vertno'].shape[0]
                # remember whether we are in the LH or RH
                this_hemi = 1
            elif ll.hemi == 'lh':
                offset = 0
                this_hemi = 0
    
            # get vertices on cortical surface inside label
            idx = np.intersect1d(ll.vertices, forward['src'][this_hemi]['vertno'])
    
            # get vertices in source space inside label
            fwd_idx = np.searchsorted(forward['src'][this_hemi]['vertno'], idx)
            thislabel_data[fwd_idx,:]=datamat[llcnt,:]
        ## project from source space to sensor space
        # multiply by leadfield
        not_noise=np.unique(np.where(thislabel_data!=0)[0])
        tobe_noise=np.asarray([thiszero for thiszero in range(thislabel_data.shape[0]) if thiszero not in not_noise])
        noise_std=np.sqrt(np.mean(RMF[125:]**2))
        noise_mat=noise_std*np.random.randn(tobe_noise.shape[0],thislabel_data[:,125:].shape[1])+np.mean(RMF[125:])
        thislabel_data[tobe_noise,125:]=noise_mat
        print ecnt
        rmat=this_std*np.random.randn(thislabel_data.shape[0],thislabel_data.shape[1])
        thislabel_datan=thislabel_data+rmat
        stc_sim_ideal[ecnt]=mne.SourceEstimate(thislabel_datan, vertices=vertices_sub, tmin=tmin, tstep=tstep, subject=subjects[subject_no], verbose=None)#thislabel_datan
        evoked_data=np.matmul(leadfield,thislabel_datan)[np.newaxis,:,:]
        if ecnt==0:
            epochs_data=evoked_data.copy()
            events=np.array([0,1,1])[np.newaxis,:]

        else:
            epochs_data=np.append(epochs_data,evoked_data,axis=0)
            this_event=np.array([ecnt+1,1,1])[np.newaxis,:]
            events=np.append(events,this_event,axis=0)
    # create evoked array
#    evoked_fwd = mne.EvokedArray(evoked_data, info=info, tmin=-0.769)
    epochs_fwd = mne.EpochsArray(epochs_data, info=epochs1.info, events=events, event_id=1, tmin=tmin)
    snr = 3.0#ediv.mean()
    lambda2 = 1.0 / snr ** 2
    method = 'MNE'
    ## apply inverse to the evoked array
    #stc_sim = apply_inverse(evoked_fwd, inverse_operator_eegmeg, lambda2,method=method, pick_ori="normal")
    stc_sim = apply_inverse_epochs(epochs_fwd, inverse_operator_eegmeg, lambda2,method=method, pick_ori="normal")
    morphmat1=mne.compute_morph_matrix(subject_from=subjects[subject_no], subject_to='fsaverage', vertices_from=stc_sim[0].vertices, vertices_to=vertices_avg, subjects_dir=data_path)
    simstc_morphed_ideal=range(len(stc_sim))
    for ccc, con1 in enumerate(stc_sim_ideal):
        print ccc		
        simstc_morphed_ideal[ccc]=con1.morph_precomputed(subject_to='fsaverage',vertices_to=vertices_avg,morph_mat=morphmat1, subject_from=subjects[subject_no])
    simstc_morphed=range(len(stc_sim))
    for ccc, con1 in enumerate(stc_sim):
        print ccc		
        simstc_morphed[ccc]=con1.morph_precomputed(subject_to='fsaverage',vertices_to=vertices_avg,morph_mat=morphmat1, subject_from=subjects[subject_no])

    # Get labels for FreeSurfer 'aparc' cortical parcellation with 34 labels/hemi
    ######################################################################################
    #load parcellatoin and do ROI connectivity
    

    label_colors = [label.color for label in labels]
    
    # Average the source estimates within each label using sign-flips to reduce
    # signal cancellations, also here we return a generator
    srcin='/imaging/rf02/TypLexMEG/fsaverage/bem/fsaverage_dist-ico-5-src.fif'
    #src_avg2=src_avg.copy()
    src_avg = mne.read_source_spaces(srcin)
    print "extracting label time course"
    label_ts = mne.extract_label_time_course(simstc_morphed, labels, src_avg, mode='mean_flip',return_generator=False)
    label_ts_ideal = mne.extract_label_time_course(simstc_morphed_ideal, labels, src_avg, mode='mean_flip',return_generator=False)
    ##################################
    
    cl=range(len(labels))
    clw=range(len(labels))
    cl_ideal=range(len(labels))
    clw_ideal=range(len(labels))
    for cntlts in range(n_epochs):
        for cntcl in range(len(labels)):        
            cl[cntcl]=simstc_morphed[cntlts].in_label(labels[cntcl]).data
            clw[cntcl]=np.argmax(np.mean(np.abs(cl[cntcl][:,125:]),axis=1))
            thislabelsign=mne.label.label_sign_flip(labels[cntcl],src_avg)[clw[cntcl]]
            label_ts[cntlts][cntcl,:]=thislabelsign*np.expand_dims(cl[cntcl][clw[cntcl],:].transpose(),0)
            ##ideal
            cl_ideal[cntcl]=simstc_morphed_ideal[cntlts].in_label(labels[cntcl]).data
            clw_ideal[cntcl]=np.argmax(np.mean(np.abs(cl_ideal[cntcl][:,125:]),axis=1))
            thislabelsign2=mne.label.label_sign_flip(labels[cntcl],src_avg)[clw_ideal[cntcl]]
            label_ts_ideal[cntlts][cntcl,:]=thislabelsign2*np.expand_dims(cl_ideal[cntcl][clw_ideal[cntcl],:].transpose(),0)
    print "extracting label time course DONE"
#############################
    
    ## correlation
#    con_corr=np.zeros((len(labels),len(labels),len(label_ts)))
#    for epochc in range(len(label_ts)):
#        con_corr[:,:,epochc]=np.corrcoef(label_ts[epochc][:,125:])
#    con_res=np.mean(con_corr,axis=2)
#    for ii1 in range(con_res.shape[0]):
#        for ii2 in range(con_res.shape[1]):
#            if ii2>=ii1:
#                con_res[ii1,ii2]=0.
#                
#    con_corr_ideal=np.zeros((len(labels),len(labels),len(label_ts_ideal)))
#    for epochc in range(len(label_ts_ideal)):
#        con_corr_ideal[:,:,epochc]=np.corrcoef(label_ts_ideal[epochc][:,125:])
#    con_res_ideal=np.mean(con_corr_ideal,axis=2)
#    for ii1 in range(con_res_ideal.shape[0]):
#        for ii2 in range(con_res_ideal.shape[1]):
#            if ii2>=ii1:
#                con_res_ideal[ii1,ii2]=0.
        
        
    # Now we are ready to compute the connectivity in the alpha band. Notice
    # from the status messages, how mne-python: 1) reads an epoch from the raw
    # file, 2) applies SSP and baseline correction, 3) computes the inverse to
    # obtain a source estimate, 4) averages the source estimate to obtain a
    # time series for each label, 5) includes the label time series in the
    # connectivity computation, and then moves to the next epoch. This
    # behaviour is because we are using generators and allows us to
    # compute connectivity in computationally efficient manner where the amount
    # of memory (RAM) needed is independent from the number of epochs.
    ##coherence
    print "computing first coh"
    fmin = 5.
    fmax = 35.
    sfreq = raw.info['sfreq']  # the sampling frequency
    con_methods = ['coh']#['pli2_unbiased']
    con, psds, freqs, times, n_epochs, n_tapers = spectral_connectivity(label_ts, method=con_methods, mode='multitaper', sfreq=sfreq, fmin=fmin,fmax=fmax, faverage=False, mt_adaptive=True, n_jobs=1)
    
    con_modified1=remove_fake_R(labels=labels, Rmat=ctf_to_ypos, con_mat=con, psd_mat=psds )
#    con_modified=np.mean(con_modified1,axis=2)
    con_modified=estimate_actual_R(labels=labels, Rmat=ctf_to_ypos, con_mat=con_modified1, psd_mat=psds )
#     con is a 3D array, get the connectivity for the first (and only) freq. band
#     for each method
    con_res = np.mean(con,axis=2)#dict()
    print "computing second coh"
    con_ideal, psds_ideal, freqs, times, n_epochs, n_tapers = spectral_connectivity(
        label_ts_ideal, method=con_methods, mode='multitaper', sfreq=sfreq, fmin=fmin,
        fmax=fmax, faverage=False, mt_adaptive=True, n_jobs=1)
    con_res_ideal = np.mean(con_ideal,axis=2)#con_ideal[:,:,0]#dict()
    print "computing 1st imcoh"
    con_methods2 = ['imcoh']#['pli2_unbiased']
    con2, psds2, freqs, times, n_epochs, n_tapers = spectral_connectivity(
        label_ts, method=con_methods2, mode='multitaper', sfreq=sfreq, fmin=fmin,
        fmax=fmax, faverage=False, mt_adaptive=True, n_jobs=1)
    con2_modified1=remove_fake_R(labels=labels, Rmat=ctf_to_ypos, con_mat=con2, psd_mat=psds )
#    con2_modified=np.mean(con2_modified1,axis=2)
    con2_modified=estimate_actual_R(labels=labels, Rmat=ctf_to_ypos, con_mat=con2_modified1, psd_mat=psds )

    
#     con is a 3D array, get the connectivity for the first (and only) freq. band
#     for each method
    con_res2 = np.mean(con2,axis=2)#con2[:,:,0]#dict()
    print "computing second imcoh"
    con_ideal2, psds_ideal2, freqs, times, n_epochs, n_tapers = spectral_connectivity(
        label_ts_ideal, method=con_methods2, mode='multitaper', sfreq=sfreq, fmin=fmin,
        fmax=fmax, faverage=False, mt_adaptive=True, n_jobs=1)
#     con is a 3D array, get the connectivity for the first (and only) freq. band
#     for each method
    con_res_ideal2 = np.mean(con_ideal2,axis=2)#con_ideal2[:,:,0]#dict()
    if cnt==0:
        con_res_final=con_res[:,:,np.newaxis].copy()
    else:
        con_res_final=np.append(con_res_final,con_res[:,:,np.newaxis],axis=2)
    if cnt==0:
        con_mod_final=con_modified[:,:,np.newaxis].copy()
    else:
        con_mod_final=np.append(con_mod_final,con_modified[:,:,np.newaxis],axis=2)
        
    if cnt==0:
        con_res_final_ideal=con_res_ideal[:,:,np.newaxis].copy()
    else:
        con_res_final_ideal=np.append(con_res_final_ideal,con_res_ideal[:,:,np.newaxis],axis=2)
    
    if cnt==0:
        con_res_final2=con_res2[:,:,np.newaxis].copy()
    else:
        con_res_final2=np.append(con_res_final2,con_res2[:,:,np.newaxis],axis=2)
    if cnt==0:
        con2_mod_final=con2_modified[:,:,np.newaxis].copy()
    else:
        con2_mod_final=np.append(con2_mod_final,con2_modified[:,:,np.newaxis],axis=2)

        
    if cnt==0:
        con_res_final_ideal2=con_res_ideal2[:,:,np.newaxis].copy()
    else:
        con_res_final_ideal2=np.append(con_res_final_ideal2,con_res_ideal2[:,:,np.newaxis],axis=2)

#    for method, c in zip(con_methods, con):
#        con_res[method] = c[:, :, 0]

#    fnamemf=data_path + meg+'simulated_sin_RMF_MTG_LO_noisefree'
#    with open(fnamemf,'wb') as f:
#        pickle.dump(stc_sim,f)                            
                            
                            
#    sim_stc=mne.simulation.simulate_stc(src_sub, morphed_labels,datamat,tmin,tstep)
    # Now, we visualize the connectivity using a circular graph layout

# First, we reorder the labels based on their location in the left hemi
label_names = [label.name for label in labels]

lh_labels = [name for name in label_names if name.endswith('lh')]

# Get the y-location of the label
label_ypos = list()
for name in lh_labels:
    idx = label_names.index(name)
    ypos = np.mean(labels[idx].pos[:, 1])
    label_ypos.append(ypos)

# Reorder the labels based on their location
lh_labels = [label for (yp, label) in sorted(zip(label_ypos, lh_labels))]

# For the right hemi
rh_labels = [label[:-2] + 'rh' for label in lh_labels]

# Save the plot order and create a circular layout
node_order = list()
node_order.extend(lh_labels[::-1])  # reverse the order
node_order.extend(rh_labels)

node_angles = circular_layout(label_names, node_order, start_pos=90,
                              group_boundaries=[0, len(label_names) / 2])

# Plot the graph using node colors from the FreeSurfer parcellation. We only
# show the 300 strongest connections.
plt.close("all")
con_res_mean=np.mean(con_res_final,axis=2)
plot_connectivity_circle(np.abs(con_res_mean), range(len(label_names)),n_lines=300,
                         node_angles=node_angles, node_colors=label_colors,#vmin=0,vmax=1,
                         title='All-to-All Connectivity '
                               'Coh')
con_mod_mean=np.mean(con_mod_final,axis=2)                             
plot_connectivity_circle(np.abs(con_mod_mean), range(len(label_names)),n_lines=300,
                         node_angles=node_angles, node_colors=label_colors,#vmin=0,vmax=1e-19,
                         title='All-to-All Connectivity '
                               'Coh')
con_res_mean_ideal=np.mean(con_res_final_ideal,axis=2)
plot_connectivity_circle(np.abs(con_res_mean_ideal), label_names, n_lines=300,#vmin=0,#vmax=1,
                         node_angles=node_angles, node_colors=label_colors,
                         title='All-to-All Connectivity '
                               'Ideal Coh')
###
con_res_mean2=np.mean(con_res_final2,axis=2)
plot_connectivity_circle(np.abs(con_res_mean2), label_names, n_lines=300,
                         node_angles=node_angles, node_colors=label_colors,#vmin=0,vmax=1,
                         title='All-to-All Connectivity '
                               'ImCoh')
con2_mod_mean=np.mean(con2_mod_final,axis=2)                             
plot_connectivity_circle(np.abs(con2_mod_mean), range(len(label_names)),n_lines=300,
                         node_angles=node_angles, node_colors=label_colors,#vmin=0,vmax=1,
                         title='All-to-All Connectivity '
                               'Coh')
con_res_mean_ideal2=np.mean(con_res_final_ideal2,axis=2)
plot_connectivity_circle(np.abs(con_res_mean_ideal2), label_names, n_lines=300,#vmin=0,vmax=1,
                         node_angles=node_angles, node_colors=label_colors,
                         title='All-to-All Connectivity '
                               'Ideal ImCoh')                               
#plt.savefig('circle.png', facecolor='black')
outname=['aparc','aparc_mod','aparc2009','aparc2009_mod','RG']
path_to=label_path+outname[which_parc]+'_coh_test.mat'#'DKA68_orig_modified_coh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot
scipy.io.savemat(path_to,{'con_res_final':con_res_final})
path_to=label_path+outname[which_parc]+'_coh_mod_test.mat'#'DKA68_orig_modified_coh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot
scipy.io.savemat(path_to,{'con_mod_final':con_mod_final})

path_to=label_path+outname[which_parc]+'_ideal_coh_test.mat'#'DKA68_ideal_modified_coh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot
scipy.io.savemat(path_to,{'con_res_final_ideal':con_res_final_ideal})

path_to=label_path+outname[which_parc]+'_imcoh_test.mat'#'DKA68_orig_modified_imcoh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot
scipy.io.savemat(path_to,{'con_res_final2':con_res_final2})
path_to=label_path+outname[which_parc]+'_imcoh_mod_test.mat'#'DKA68_orig_modified_imcoh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot
scipy.io.savemat(path_to,{'con2_mod_final':con2_mod_final})

path_to=label_path+outname[which_parc]+'_ideal_imcoh_test.mat'#'DKA68_ideal_modified_imcoh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot
scipy.io.savemat(path_to,{'con_res_final_ideal2':con_res_final_ideal2})
plt.show()
        