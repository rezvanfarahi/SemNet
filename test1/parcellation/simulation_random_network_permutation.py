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
#import sklearn
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

#subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
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
which_parc=0 #names below
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
n_permutation=100
n_subj=len(ll)
coh_ideal_max=np.zeros((n_subj,n_permutation))
imcoh_ideal_max=np.zeros((n_subj,n_permutation))
for cnt, meg in enumerate(ll):
    ii=ii+1
    print cnt
    subjects_dir = data_path 
    fname_fwd = data_path + meg + 'ico5_forward_5-3L-EMEG-fwd.fif'
    fname_inv = inv_path + meg + 'InvOp_ico5newreg_fft_1_48_clean_ica_EMEG-inv.fif'
    vertices_avg = [np.arange(10242), np.arange(10242)]
    subject_no=subject_inds[0]
    raw_fname = data_path + meg + 'semloc_ssstf_fft_1_48_clean_ica_raw.fif' #clean_ssp_
    event_fname = event_path + meg + 'semloc_raw_ssstnew.txt'
    subjects_dir = data_path 
    
    
    raw = io.Raw(raw_fname)
    events1 = mne.read_events(event_fname)
    stim_delay = 0.034 # delay in s
    events1[:,0]= events1[:,0]+np.round( raw.info['sfreq']*stim_delay )
    
    # Read epochs for all channels, removing a bad one
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
    
    for rndcnt in range(n_permutation):
        print rndcnt
        n_epochs=50
        LO_phase=np.pi*(2*np.random.random(n_epochs,)-1)
        RMF_phase=np.pi*(2*np.random.random(n_epochs,)-1)
        stc_sim_ideal=range(n_epochs)
        for ecnt in range(n_epochs):
            t=np.linspace(-6*np.pi,6*np.pi,3020/5)
            RMF=np.concatenate((np.zeros(125,),1e-9*np.sin(2*t+RMF_phase[ecnt])),axis=0)#+RMF_phase[ecnt]
            this_std=np.sqrt(np.mean(RMF[125:]**2)/9.)
            tmin=-0.124
            tstep=0.001
            lvert=src_sub[0]['vertno']
            rvert=src_sub[1]['vertno']
            vertices_sub = [lvert, rvert]
            thislabel_data=np.zeros((len(lvert)+len(rvert),len(RMF)))
            tobe_noise=np.arange(thislabel_data.shape[0])#np.asarray([thiszero for thiszero in range(thislabel_data.shape[0]) if thiszero not in not_noise])
            noise_std=np.sqrt(np.mean(RMF[125:]**2))
            noise_mat=noise_std*np.random.randn(tobe_noise.shape[0],thislabel_data[:,125:].shape[1])+np.mean(RMF[125:])
            thislabel_data[tobe_noise,125:]=noise_mat
            #            print ecnt
            rmat=this_std*np.random.randn(thislabel_data.shape[0],thislabel_data.shape[1])
            thislabel_datan=thislabel_data+rmat
            stc_sim_ideal[ecnt]=mne.SourceEstimate(thislabel_datan, vertices=vertices_sub, tmin=tmin, tstep=tstep, subject=subjects[subject_no], verbose=None)#thislabel_datan
        snr = 3.0#ediv.mean()
        lambda2 = 1.0 / snr ** 2
        method = 'MNE'
        ## apply inverse to the evoked array
        #stc_sim = apply_inverse(evoked_fwd, inverse_operator_eegmeg, lambda2,method=method, pick_ori="normal")
#        sub_vertices=[inverse_operator_eegmeg['src'][0]['vertno'],inverse_operator_eegmeg['src'][1]['vertno']]
        morphmat1=mne.compute_morph_matrix(subject_from=subjects[subject_no], subject_to='fsaverage', vertices_from=vertices_sub, vertices_to=vertices_avg, subjects_dir=data_path)
        simstc_morphed_ideal=range(n_epochs)
        for ccc, con1 in enumerate(stc_sim_ideal):
#            print ccc		
            simstc_morphed_ideal[ccc]=con1.morph_precomputed(subject_to='fsaverage',vertices_to=vertices_avg,morph_mat=morphmat1, subject_from=subjects[subject_no])
    
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
        label_ts_ideal = mne.extract_label_time_course(simstc_morphed_ideal, labels, src_avg, mode='mean_flip',return_generator=False)
        ##################################
        
        cl_ideal=range(len(labels))
        clw_ideal=range(len(labels))
        for cntlts in range(n_epochs):
            for cntcl in range(len(labels)):        
                ##ideal
                cl_ideal[cntcl]=simstc_morphed_ideal[cntlts].in_label(labels[cntcl]).data
                clw_ideal[cntcl]=np.argmax(np.mean(np.abs(cl_ideal[cntcl][:,125:]),axis=1))
                thislabelsign2=mne.label.label_sign_flip(labels[cntcl],src_avg)[clw_ideal[cntcl]]
                label_ts_ideal[cntlts][cntcl,:]=thislabelsign2*np.expand_dims(cl_ideal[cntcl][clw_ideal[cntcl],:].transpose(),0)
        print "extracting label time course DONE"
        #############################
            
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
        fmin = 5.
        fmax = 35.
        sfreq = raw.info['sfreq']  # the sampling frequency
        con_methods = ['coh']#['pli2_unbiased']
        print "computing ideal coh"
        con_ideal, psds_ideal, freqs, times, n_epochs, n_tapers = spectral_connectivity(
            label_ts_ideal, method=con_methods, mode='multitaper', sfreq=sfreq, fmin=fmin,
            fmax=fmax, faverage=True, mt_adaptive=True, n_jobs=4)
        con_res_ideal = con_ideal[:,:,0]#np.mean(con_ideal,axis=2)#con_ideal[:,:,0]#dict()
        coh_ideal_max[cnt,rndcnt]=np.max(np.abs(con_res_ideal))
        print "computing ideal imcoh"
        con_methods2 = ['imcoh']
        con_ideal2, psds_ideal2, freqs, times, n_epochs, n_tapers = spectral_connectivity(
            label_ts_ideal, method=con_methods2, mode='multitaper', sfreq=sfreq, fmin=fmin,
            fmax=fmax, faverage=True, mt_adaptive=True, n_jobs=4)
    #     con is a 3D array, get the connectivity for the first (and only) freq. band
    #     for each method
        con_res_ideal2 = con_ideal2[:,:,0]# np.mean(con_ideal2,axis=2)#con_ideal2[:,:,0]#dict()
        imcoh_ideal_max[cnt,rndcnt]=np.max(np.abs(con_res_ideal2))
        
    outname=['aparc_','aparc_mod_','aparc2009_','aparc2009_mod_','RG_']
    path_to=label_path+outname[which_parc]+subjects[subject_inds[0]]+'_coh_ideal_max.mat'#'DKA68_ideal_modified_coh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot
    scipy.io.savemat(path_to,{'coh_ideal_max':coh_ideal_max})
    
    path_to=label_path+outname[which_parc]+subjects[subject_inds[0]]+'_imcoh_ideal_max.mat'#'DKA68_orig_modified_imcoh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot
    scipy.io.savemat(path_to,{'imcoh_ideal_max':imcoh_ideal_max})
