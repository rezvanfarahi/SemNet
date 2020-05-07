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

# Authors: Rezvan Farahibozorg 09.08.2017
#
# License: BSD (3-clause)



# Author: Martin Luessi <mluessi@nmr.mgh.harvard.edu>
#
# License: BSD (3-clause)

print(__doc__)
print "this script is to reconstruct simulated E/MEG connectomes in order to compare the performance of coherence and multivariate coherence"
import sys

#sys.path.append('/imaging/local/software/python_packages/nibabel/1.3.0')
#sys.path.append('/imaging/local/software/python_packages/pysurfer/v0.4')
#sys.path.append('/imaging/rf02/scikit-learn-0.15.0')
#
##sys.path.append('/imaging/local/software/python_packages/scikit-learn/v0.14.1')
#sys.path.append('/imaging/local/software/mne_python/la_v0.9')
#sys.path.insert(1,'/imaging/rf02/mne_python_11')
#sys.path.insert(1,'/imaging/local/software/mne_python/la_v0.9full')
sys.path.insert(1,'/imaging/local/software/mne_python/latest_v0.15')


sys.path.insert(1,'/imaging/rf02/scikit-learn-0.16.1')
#sys.path.insert(1,'/home/rf02/.local/lib/python2.7/site-packages')


# for qsub
# (add ! here if needed) 
#sys.path.insert(1,'/imaging/local/software/EPD/latest/x86_64/bin/python')
#sys.path.insert(1,'/imaging/local/software/anaconda/2.4.1/2/lib/python2.7/site-packages')
#sys.path.insert(1,'/imaging/local/software/anaconda/latest/x86_64/pkgs/python-2.7.13-0')
#sys.path.insert(1,'/imaging/local/software/anaconda/latest/x86_64')
#sys.path.insert(1,'/imaging/local/software/anaconda/latest/3/lib/python3.5/site-packages')
sys.path.insert(1,'/imaging/local/software/anaconda/2.4.1/2/lib/python2.7')
#sys.path.append('/imaging/rf02/TypLexMEG/mne_python_v8/')
#sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
###
#import sklearn
#from scipy.stats.mstats import zscore
import os
import numpy as np
#reload(np)
print np.__path__
print np.__version__
import mne
reload(mne)
#print mne.__path__
from mne import io
#from mne.io import Raw
#import operator
from mne.minimum_norm import (read_inverse_operator, apply_inverse_epochs)
import scipy.io
from mne.label import _n_colors as nc
from mne.label import _read_annot as mnera
from mne.label import _write_annot as mnewa
from mne.connectivity import spectral_connectivity
from mne.epochs import equalize_epoch_counts

import matplotlib.pyplot as plt
import copy
import pickle
print mne.__path__
import time
import math


data_path = '/imaging/rf02/Semnet/' # root directory for your MEG data
os.chdir(data_path)
event_path = '/imaging/rf02/Semnet/'#'/group/erp/data/olaf.hauk/MEG/TypLexM/data/MF22_new/'
#label_path = '/group/erp/data/olaf.hauk/MEG/TypLexM/labels/createdlabels_SL/'
inv_path = '/imaging/rf02/Semnet/'
subjects_dir = data_path 
print sys.argv
subject_inds = []#
for ss in sys.argv[1:]:
    subject_inds.append( int( ss ) )

subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]
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
for ss in subject_inds:
 ll.append(list_all[ss])


print "ll:"
print ll

label_path='/imaging/rf02/TypLexMEG/parcellation/labels_and_CTFs/'
results_path='/imaging/rf02/parcellation/leakage_pattern/'
if not os.path.exists(results_path):
    os.makedirs(results_path)

##grow ROIs
srcin='/imaging/rf02/TypLexMEG/fsaverage/bem/fsaverage_dist-ico-5-src.fif'
#src_avg2=src_avg.copy()
src_avg = mne.read_source_spaces(srcin)
subject_avg='fsaverage'
#nseeds=20


#for lbl in labellist_m:
#    lbl.values.fill(1.0)
##
#which_parc=1 #names below
parc_name=['aparc',
'fsaverage_finalpostparc_smooth_stddiff_split_merge_halfmax_86_badremoved_74',
'aparc.a2009s',
'fsaverage_finalpostparc2009_smooth_stddiff_split_merge_halfmax_98_badremoved_74_cmatch',
'fsaverage_RG_finalpostparc_smooth_stddiff_split_merge_halfmax_82_badremoved_70_cmatch'
]##'fsaverage_finalpostparc2009_smooth_stddiff_split_merge_halfmax_98_badremoved_74_cmatch'#'fsaverage_RG_finalpostparc_smooth_stddiff_split_merge_halfmax_82_badremoved_70_cmatch'
#
Rmat_names=['ctf_ypos_aparc.mat',
            'ctf_ypos_finalpostparc_smooth_stddiff_split_merge_halfmax_86_badremoved_74.mat',
            'ctf_ypos_anat2009.mat',
            'ctf_ypos_finalpostparc2009_smooth_stddiff_split_merge_halfmax_98_badremoved_74_cmatch.mat',
            'ctf_ypos_RG_finalpostparc_smooth_stddiff_split_merge_halfmax_82_badremoved_70_cmatch.mat'
            ]

ypos_names=['ypos_aparc.mat',
            'ypos_aparcmod_badremoved_74.mat',
            'ypos_aparc2009.mat',
            'ypos_aparc2009mod_badremoved_74.mat',
            'ypos_RG_badremoved_70.mat'
            ]


parc_len=np.zeros((len(parc_name),))
sldminmean=np.zeros((len(parc_name),))
for which_parc in range(len(parc_name)):
    labels_len = mne.read_labels_from_annot('fsaverage', parc=parc_name[which_parc], subjects_dir=subjects_dir)
    if which_parc==0:
        labels_len.remove(labels_len[-1])
    if which_parc==2:
        labels_len.remove(labels_len[-1])
        labels_len.remove(labels_len[-1])
    parc_len[which_parc]=len(labels_len)
               

outparc=['aparc','aparc_mod','aparc2009','aparc2009_mod','RG']
           
#labels_pre.remove(labels_pre[-1])
#    labels_pre.remove(labels_pre[-1])
for cnt, meg in enumerate(ll):
    print cnt
    subject_no=subject_inds[cnt]
    subjects_dir = data_path 
    fname_fwd = data_path + meg + 'ico5_forward_5-3L-EMEG-fwd.fif'
    forward = mne.read_forward_solution(fname_fwd,surf_ori=True, force_fixed=True)#, surf_ori=True)
    src_sub=forward['src']
    leadfield=forward['sol']['data']
    lvert=src_sub[0]['vertno']
    rvert=src_sub[1]['vertno']
    vertices_sub = [lvert, rvert]
#    fname_inv = inv_path + meg + 'InvOp_ico5newreg_filtnew_pnt1_48_clean_ica_EMEG-inv.fif'#'InvOp_ico5diagreg_fft_1_48_clean_ica_EMEG-inv.fif'
    vertices_avg = [np.arange(10242), np.arange(10242)]
    lvertavg=src_avg[0]['vertno']
    rvertavg=src_avg[1]['vertno']
    
    srcin='/imaging/rf02/TypLexMEG/fsaverage/bem/fsaverage_dist-ico-5-src.fif'
    src = mne.read_source_spaces(srcin)
    subjects_dir = data_path 
    print "load raw"
    raw_fname = data_path + meg+ 'SemDec_blocks_tsss_filtnew_pnt1_48_ica_raw.fif' #clean_ssp_
    event_fname = data_path + meg + 'SemDec_blocks_tsss_filtnew_pnt1_48_ica_raw-eve.fif'
    fname_stcforverts = data_path + meg + 'firstMorphed_ico_SemDec_pnt1_30ica_Visual_Source_Evoked_m300_600'#_'+b#'firstMorphed_ico_SemDec_ica_'+event_names[event_no]+'_Source_Evoked_mnrsmpl_50_550_100ms_0ov'#'_Source_Evoked_m300_600'#'_Source_Evoked_mnrsmpl_50_550_100ms_0ov'
    stc_forverts = mne.read_source_estimate(fname_stcforverts)
    subjects_dir = data_path 
    
    tmin = -0.3
    tmax = 0.6  # Use a lower tmax to reduce multiple comparisons
    #   Setup for reading the raw data
    raw = io.Raw(raw_fname, preload=True)
    events = mne.read_events(event_fname)
    for evcnte in range(events.shape[0]-1):
        if events[evcnte,2] in np.array([2]):
            if events[evcnte+1,2] in np.array([79,90]):
                events[evcnte,2]=230
        if events[evcnte,2] in np.array([1]):
            if events[evcnte+1,2] in np.array([71]):
                events[evcnte,2]=230
    #    stim_delay = 0.034 # delay in s
    #    events[:,0] = events[:,0]+np.round( raw.info['sfreq']*stim_delay )
    ###############################################################################
    # Read epochs for all channels, removing a bad one
    picks = mne.pick_types(raw.info, eeg=True, meg=True, eog=False, exclude='bads')
    print "Epoching"
    if subject_no in np.array([2,3,5,18]):#np.array([4,8,9]):
        reject = dict(eeg=150e-6, grad=200e-12, mag=4e-12)
    else:
        reject = dict(eeg=120e-6, grad=200e-12, mag=4e-12)#, eog=150e-6)
    event_ids = {'visual': 1, 'hear': 2, 'hand': 3, 'pwordc': 6}#'neutral': 4, 'emotional': 5,
    epochs = mne.Epochs(raw, events, event_ids, tmin, tmax, picks=picks, proj=True, baseline=(None, 0), reject=reject,preload=True)
    epochs.drop_bad()
    print len(epochs)/600.
    epochs_list=range(4)
    for eid,event_id in enumerate(list(np.array([1,2,3,6]))):
#                event_id = eid+1  
        epochs_list[eid]=epochs[epochs.events[:,2]==event_id]
    equalize_epoch_counts(epochs_list,method='mintime')
    print "Transform to source space"

    fname_inv = data_path + meg + 'InvOp_ico5newreg_filtnew_pnt1_48_clean_ica_EMEG-inv.fif'

    method = "MNE"  # use dSPM method (could also be MNE or sLORETA)
    inverse_operator = read_inverse_operator(fname_inv)

    subject_from = subjects[subject_no]
    subject_to = 'fsaverage'
    vertices_to = [np.arange(10242), np.arange(10242)]
    event_names = ['Visual', 'Hear', 'Hand', 'Pwordc']#'Neutral', 'Emotional',
    
    for which_parc in np.array([1]):#range(len(parc_name)):#np.array([2,3]):#
        print outparc[which_parc]
        labels_pre = mne.read_labels_from_annot('fsaverage', parc=parc_name[which_parc], subjects_dir=subjects_dir)
        if which_parc==0:
            labels_pre.remove(labels_pre[-1])
        if which_parc==2:
            labels_pre.remove(labels_pre[-1])
            labels_pre.remove(labels_pre[-1])
        #labels_pre.remove(labels_pre[-1])
        #    labels_pre.remove(labels_pre[-1])
        
        path_Rmat=results_path+Rmat_names[which_parc]
        ltc=scipy.io.loadmat(path_Rmat,mat_dtype=True); ctf_to_ypos=np.squeeze(ltc['ctf_to_ypos'])            
        
        path_ypos=results_path+ypos_names[which_parc]
        ltc=scipy.io.loadmat(path_ypos,mat_dtype=True); ypos_sort=np.squeeze(ltc['ypos_sort'])
        labellist=[labels_pre[jj] for jj in ypos_sort]
        for lbl in labellist:
            lbl.values.fill(1.0)
        labels=list()
        
        ctfstc_path=data_path+'/stc_ctf_meaneve_'+outparc[which_parc]
        stc_subsm = mne.read_source_estimate(ctfstc_path)
        Matx=stc_subsm.data[:,:-1].copy()
        print "pick label verts"
        for this_labelc, this_label in enumerate(labellist): 
            vertices_test = [np.arange(10242), np.arange(10242)]

            label_origdata=stc_forverts.in_label(this_label).data
            if this_label.hemi=='lh':
                label_verts=stc_forverts.in_label(this_label).lh_vertno
                vert50=np.sort(Matx[label_verts,this_labelc])[-1]
                winner_verts=label_verts[np.where(Matx[label_verts,this_labelc]>=vert50)[0]]#label_verts[np.where(Matx[label_verts,this_labelc]>=np.max(Matx[label_verts,this_labelc])/2.)[0]]#np.where(Matx[:,this_labelc]>=np.max(Matx[:,this_labelc])/2.)[0]#

            else:
                label_verts=stc_forverts.in_label(this_label).rh_vertno+len(vertices_test[0])
                vert50=np.sort(Matx[label_verts,this_labelc])[-1]
                winner_verts=label_verts[np.where(Matx[label_verts,this_labelc]>=vert50)[0]]#label_verts[np.where(Matx[label_verts,this_labelc]>=np.max(Matx[label_verts,this_labelc])/2.)[0]]#np.where(Matx[:,this_labelc]>=np.max(Matx[:,this_labelc])/2.)[0]#

            #            np.intersect1d(label_verts,vertices_to[0])
            #            #            print winner_verts
#            thisscale=np.hstack((mne.label.label_sign_flip(thisstc_label[0],src_avg),mne.label.label_sign_flip(thisstc_label[1],src_avg)))[winner_verts]#mne.label.label_sign_flip(this_label,src_avg)[winner_verts]#1#Matx[winner_verts,this_labelc]/np.max(Matx[winner_verts,this_labelc])
#            label_predata=this_stc.data[winner_verts,:]
#            #            label_data1=np.tile(thisscale,[ntimes,1]).transpose()*label_predata
#            label_data1=thisscale[:,np.newaxis]*label_predata
#            active_ones=np.where(np.mean(np.abs(label_data1[:,sonset:]),1)>=np.min(np.min(np.abs(label_data1[:,sonset:]),1)))[0]
#            print len(active_ones)
#            label_data=label_data1[active_ones,:]
#            ini_tc=np.mean(label_data[:,:], axis=0)
            thisstc_testdata=np.zeros((stc_forverts.data.shape[0],1))
            thisstc_testdata[winner_verts,0]=1
            thisstc_tolabel_test = mne.SourceEstimate(thisstc_testdata, vertices=vertices_test,tmin=1e-3 , tstep=1e-3 , subject='fsaverage')
            if this_label.hemi=='lh':
                thisstc_labeltest=mne.stc_to_label(thisstc_tolabel_test, src_avg,smooth=False)[0]
            else:
                thisstc_labeltest=mne.stc_to_label(thisstc_tolabel_test, src_avg,smooth=False)[1]

            thislabel_path=data_path+'ROI_analysis/parcel_connectivity/'+outparc[which_parc]+'ctflabs/'
            if not os.path.exists(thislabel_path):
                os.makedirs(thislabel_path)
            testlabel_path=thislabel_path+'subject'+str(subject_no)+'_'+this_label.name
            thisstc_labeltest.save(testlabel_path)
            labels.append(thisstc_labeltest)
            #            conabs_mat[this_labelc,:,event_no]=ini_tc.copy()#(ini_tc-np.min(ini_tc))/(np.max(ini_tc)-np.min(ini_tc))zscore(ini_tc)#


        for evcnt, this_epoch in enumerate(epochs_list):
            print event_names[evcnt]; print "apply inverse epochs"
            print len(this_epoch)/150.
#            this_epoch = copy.deepcopy(epochs_list[evcnt])

            snr = 1.0#np.mean(ediv)
            lambda2 = 1.0 / snr ** 2
            this_condition = apply_inverse_epochs(this_epoch, inverse_operator, lambda2,method=method, 
            pick_ori="normal")
            morphmat_sub2avg=mne.compute_morph_matrix(subject_from=subjects[subject_no], 
            subject_to='fsaverage', vertices_from=vertices_sub, vertices_to=vertices_avg, 
            subjects_dir=data_path)
            condition_morphed=[con1.copy().morph_precomputed(subject_to='fsaverage',vertices_to=vertices_avg,morph_mat=morphmat_sub2avg, subject_from=subjects[subject_no]) for con1 in this_condition]
            print "extracting label time course"
            label_ts = mne.extract_label_time_course(condition_morphed, labels, src_avg, mode='mean_flip',return_generator=False)
            print "computing connectivity"
            fmin = (4.,8.,12.,30.)
            fmax = (8.,12.,30.,45.)
            sfreq = raw.info['sfreq']  # the sampling frequency
            con_methods = ['coh']#['pli2_unbiased']
            connmat=np.zeros((len(labels),len(labels),len(fmin),3))#np x np x nf x nt
            tcnt=-1
            for tminc, tmaxc in zip([0.05,0.15,0.25],[0.25,0.35,0.45]):
                tcnt=tcnt+1
                con, freqs, times, n_epochs, n_tapers = spectral_connectivity(
                label_ts, method=con_methods, mode='multitaper', sfreq=sfreq, fmin=fmin,fmax=fmax, 
                tmin=tminc,tmax=tmaxc,faverage=True, mt_adaptive=True, n_jobs=1)                
                connmat[:,:,:,tcnt] = np.squeeze(con.copy())#np.mean(con,axis=2)#dict()
            outparc=['aparc','aparc_mod','aparc2009','aparc2009_mod','RG']
            out_path=data_path+meg+'pacel_connmat_maxctf_'+event_names[evcnt]+'_'+outparc[which_parc]+'.mat'
            scipy.io.savemat(out_path,{'connmat':connmat})


####  HERE!!!
#
#        for event_no, event_name in enumerate(event_names):
#            for b in bands:
#                fname_in = data_path + meg + 'firstMorphed_SemDec_filtnew_pnt1_48ica_'+event_names[event_no]+'_Source_signedEvoked_m300_600_'+b#'firstMorphed_ico_SemDec_pnt1_30ica_'+event_names[event_no]+'_Source_Evoked_mnrsmpl_50_550_100ms_0ov'#'_Source_Evoked_m300_600'#'_Source_Evoked_mnrsmpl_50_550_100ms_0ov'
#                this_stc_normal = mne.read_source_estimate(fname_in)
#    
#    
#                label_ts_normal = mne.extract_label_time_course(this_stc_normal, labellist, src, mode='mean_flip',return_generator=False)
#                stc_sim = apply_inverse_epochs(epochs_fwd, inverse_operator_eegmeg, lambda2,method=method, pick_ori="normal")
#
#
#    print "epochs"
#    picks = mne.pick_types(raw.info, eeg=True, meg=True, eog=False, exclude=[])
#    epochs1 = mne.Epochs(raw, events1, tmin=-0.5, tmax=0.7,event_id=1, picks=picks, baseline=(None, 0), proj=True, reject=None, preload=True)
#    
#    
#    # read inverse operators
#    inverse_operator_eegmeg = read_inverse_operator(fname_inv)
#    info = copy.deepcopy(inverse_operator_eegmeg['info'])
#    info['sfreq'] = 1000.  # necessary
#    info['projs'] = inverse_operator_eegmeg['projs']
#    ## simulate activity in labels
#                
#    src_sub=forward['src']
#    leadfield=forward['sol']['data']
#    lvert=src_sub[0]['vertno']
#    rvert=src_sub[1]['vertno']
#    vertices_sub = [lvert, rvert]
#    
#    
#        labellist_s=copy.deepcopy(labellist)
#        labov_mat=np.zeros((len(labellist_s),len(labels)))
#        labels_mod=copy.deepcopy(labels)
#        lpcnt=-1
#        for lscnt,seedlab in enumerate(labellist_s):
#            for lpcnt,parclab in enumerate(labels):
#                if seedlab.name[-2:]==parclab.name[-2:]:
#                    labov_mat[lscnt,lpcnt]= len(np.intersect1d(seedlab.vertices, parclab.vertices))
##                    else:
##                        labov_mat[lscnt,lpcnt]=0                  
#        iarpre=np.argmax(labov_mat,1)
#        im=np.max(labov_mat,1)
#        iar=iarpre.copy()
#        iarpos=iarpre[im>0]
#        iar[im==0]=1000
#        misscatch=iarpre.copy()
#        misscatch[im==0]=0
#        misscatch[im>0]=1
#        print iar
#        missed_src=len(iarpre)-len(iar)
#        actsrc_miss=iar.copy()
#        labov_mat_final=labov_mat.copy()
#        labov_mat_final[labov_mat_final>0]=1
#        ov_sum=np.sum(labov_mat_final,1)
#        ov_sum[ov_sum<1]=1
#        howmany_actsrc=np.sum(ov_sum)
#        actsrc_final=np.where(labov_mat_final==1)[1]
##                morphed_labels=[labellist_s[jj].morph(subject_from='fsaverage', smooth=1, subject_to=subjects[subject_no], subjects_dir=data_path) for jj in range(len(labellist_s))]
#        
#
##                fname_split=parc_name[which_parc]+'_extsims'#+str(perm_inds[0])
##                
##                morphed_labelsall = mne.read_labels_from_annot(subject=subjects[subject_no], parc=subjects[subject_no]+'_'+fname_split, subjects_dir=data_path)
##
##                print "hi"
##                for lblcnt,labelcheck in enumerate(morphed_labelsall):
##                    if labelcheck.hemi=='lh':
##                        lblcs=np.intersect1d(src_sub[0]['vertno'],labelcheck.vertices)
##                    else:
##                        lblcs=np.intersect1d(src_sub[1]['vertno'],labelcheck.vertices)
##                    if len(lblcs)<1:
##                        print labelcheck
##                        labellist_mc=copy.deepcopy(labels)
##                        morphed_labelsall[lblcnt]=labellist_mc[lblcnt].morph(subject_from='fsaverage', smooth=5, subject_to=subjects[subject_no], subjects_dir=data_path)
#
##                labelsp=copy.deepcopy(morphed_labelsall)
##                pick_seeds=[labelsp[pscnt] for pscnt in iarpos]#fsaverage space
##    #            morphed_labels=copy.deepcopy(pick_seeds)
##                for thisps in pick_seeds:
##                    thisps.values.fill(1.0)
##                pick_seeds_cmas=np.asarray([llcm.center_of_mass(subject=subjects[subject_no],restrict_vertices=True, subjects_dir=data_path) for llcm in pick_seeds])
##                pick_seeds_heml=list()
##                for pscnt,llcm in enumerate(pick_seeds):
##                    if llcm.hemi=='lh':
##                        pick_seeds_heml.append(0)
##                    else:
##                        pick_seeds_heml.append(1)
##                pick_seeds_hem=np.asarray(pick_seeds_heml)
##                final_ext=5
##                morphed_labels=mne.grow_labels(subject=subjects[subject_no],seeds=pick_seeds_cmas,extents=final_ext,hemis=pick_seeds_hem,subjects_dir=data_path,overlap=False)
#
##            morphed_labels=[picksds[jj].morph(subject_from='fsaverage', smooth=5, subject_to=subjects[subject_no], subjects_dir=data_path) for jj in range(len(picksds))]
##                morphed_labels=copy.deepcopy(labellist)
#        morphed_labelsall=copy.deepcopy(labels)
##                labov_mat_final=np.zeros((len(morphed_labels),len(morphed_labelsall)))
##                lpcnt1=-1
##                for lscnt,seedlab in enumerate(morphed_labels):
##                    for lpcnt1,parclab in enumerate(morphed_labelsall):
##                        if seedlab.name[-2:]==parclab.name[-2:]:
##                            labov_mat_final[lscnt,lpcnt1]= len(np.intersect1d(seedlab.vertices, parclab.vertices))
##                labov_mat_final[labov_mat_final>0]=1
##                ov_sum=np.sum(labov_mat_final,1)
##                ov_sum[ov_sum<1]=1
##                howmany_actsrc=np.sum(ov_sum)
##                actsrc_final=np.where(labov_mat_final==1)[1]
#        
#        #labov_mat[labov_mat>0]=1
#        idealconn=np.eye(len(morphed_labelsall))
#        amatt=Amat.copy()
#        amatt=amatt[misscatch>0,:]
#        amatt=amatt[:,misscatch>0]
#        if len(iarpos)>1:
#            for ccc1,iarc1 in enumerate(iarpos):
#                for ccc2,iarc2 in enumerate(iarpos):
#                    if amatt[ccc1,ccc2]>0:
#                        idealconn[iarc1,iarc2]=1
#        for lovc in range(labov_mat_final.shape[0]):
#            this_row=np.where(labov_mat_final[lovc,:]>0)[0]
#            if len(this_row)>1:
#                for ccc1,iarc1 in enumerate(this_row):
#                    for ccc2,iarc2 in enumerate(this_row):
#                        idealconn[iarc1,iarc2]=1
#        for lovc in range(labov_mat_final.shape[0]):
#            this_row=np.where(labov_mat_final[lovc,:]>0)[0]
#            for this_roi1 in this_row:
#                for this_roi2 in this_row:
#                    for this_parcel in range(idealconn.shape[0]):
#                        if idealconn[this_roi1,this_parcel]==1:# or idealconn[this_parcel,this_roi1]>0:
#                            idealconn[this_roi2,this_parcel]=1
#                            idealconn[this_parcel,this_roi2]=1
#        for lovc in range(labov_mat_final.shape[0]):
#            this_row=np.where(labov_mat_final[lovc,:]>0)[0]
#            this_roi1=iar[lovc].copy()
#            for this_roi2 in this_row:
#                for ii1 in range(idealconn.shape[1]):
#                    if idealconn[this_roi1,ii1]==1:
#                        idealconn[this_roi2,ii1]=1
####NOT SURE IF USED BELOW OR ABOVE IN THE END, guessing above- DOUBLE check...                               
##                idealconn2=np.eye(len(morphed_labelsall))
##                amatt=Amat.copy()
##                amatt=amatt[misscatch>0,:]
##                amatt=amatt[:,misscatch>0]
##                if len(iarpos)>1:
##                    for ccc1,iarc1 in enumerate(iarpos):
##                        for ccc2,iarc2 in enumerate(iarpos):
##                            if amatt[ccc1,ccc2]>0:
##                                idealconn2[iarc1,iarc2]=1
##                                idealconn2[iarc2,iarc1]=1
##                for lovc in range(labov_mat_final.shape[0]):
##                    this_row=np.where(labov_mat_final[lovc,:]>0)[0]
##                    if len(this_row)>1:
##                        for ccc1,iarc1 in enumerate(this_row):
##                            for ccc2,iarc2 in enumerate(this_row):
##                                idealconn2[iarc1,iarc2]=1
##                for lovc in range(labov_mat_final.shape[0]):
##                    this_row=np.where(labov_mat_final[lovc,:]>0)[0]
##                    for this_roi1 in this_row:
##                        for this_roi2 in this_row:
##                            for this_parcel in range(idealconn2.shape[0]):
##                                if idealconn2[this_roi1,this_parcel]==1:# or idealconn2[this_parcel,this_roi1]>0:
##                                    idealconn2[this_roi2,this_parcel]=1
##                                    idealconn2[this_parcel,this_roi2]=1
##                for lovc in range(labov_mat_final.shape[0]):
##                    this_row=np.where(labov_mat_final[lovc,:]>0)[0]
##                    this_roi1=iar[lovc].copy()
##                    for this_roi2 in this_row:
##                        for ii1 in range(idealconn2.shape[1]):
##                            if idealconn2[this_roi1,ii1]==1:
##                                idealconn2[this_roi2,ii1]=1
##                print "terminate here"
#                        
#            
#            
#        extramisses=0
#        for seed_pairs in comseed:
##                    idealconn[seed_pairs[0],seed_pairs[1]]=1
##                    idealconn[seed_pairs[1],seed_pairs[0]]=1
#            if seed_pairs[0] in np.where(misscatch==0)[0] or seed_pairs[1] in np.where(misscatch==0)[0]:
#                extramisses=extramisses+1
#        missedcons=np.tril(Amat)-np.eye(Amat.shape[0])
#        #            missedcons=missedcons[misscatch==0,:]
#        mct=missedcons[misscatch==0,:]
#        mct=mct[:,misscatch==0]
#        missedcons=np.sum(missedcons[misscatch==0,:])+np.sum(missedcons[:,misscatch==0])-np.sum(mct)+extramisses
#        if cnt==0:
#            actsrc_grandmat.append(actsrc_miss)
#            actsrcfinal_grandmat.append(actsrc_final)
#            missedcons_grandmat.append(missedcons)
#            actsrc_ovs.append(ov_sum)
#            misscatch_grandmat.append(misscatch)
##            labov_mat[labov_mat>0]=1
##            ov_sum=np.sum(labov_mat,1)
##            ov_sum[ov_sum<1]=1
##            howmany_actsrc=np.sum(ov_sum)
##            labellist_m=copy.deepcopy(pick_seeds)
##[labellist_m[jj].morph(subject_from='fsaverage', smooth=5, subject_to=subjects[subject_no], subjects_dir=data_path) for jj in range(len(labellist_m))]
#        parc_datamat_grand=datamat_grand.copy()#datamat_grand[np.where(im>0)[0],:,:].copy()
##                lowsnr_parc_datamat_grand=lowsnr_datamat_grand.copy()#lowsnr_datamat_grand[np.where(im>0)[0],:,:].copy()
#        label_colors = [label.color for label in labels]
#        stc_sim_ideal=range(n_epochs)
#        stc_sim_realideal=range(n_epochs)
##                lowsnr_stc_sim=range(n_epochs)
##                lowsnr_stc_sim_ideal=range(n_epochs)
#        for ecnt in range(n_epochs):  
#            datamat=parc_datamat_grand[:,:,ecnt].copy()
##                    lowsnr_datamat=lowsnr_parc_datamat_grand[:,:,ecnt].copy()
#            this_std=np.sqrt(np.mean(datamat[0,blen:]**2)/(thissnr**2))#SNR=3 
##                    lowsnr_this_std=np.sqrt(np.mean(lowsnr_datamat[0,blen:]**2)/(lowsnr**2))
#            tmin=-0.124
#            tstep=0.001
#            
#            thislabel_data=np.zeros((len(lvertavg)+len(rvertavg),datamat.shape[1]))
##                    allnoise_thislabel_data=np.zeros((len(lvert)+len(rvert),datamat.shape[1]))
##                    lowsnr_thislabel_data=np.zeros((len(lvertavg)+len(rvertavg),lowsnr_datamat.shape[1]))
#       
#            for llcnt,ll1 in enumerate(labellist):
#                if ll1.hemi == 'rh':
#                    # for RH labels, add number of LH vertices
#                    offset = len(lvertavg)#src_avg[0]['vertno'].shape[0]
#                    # remember whether we are in the LH or RH
#                    this_hemi = 1
#                elif ll1.hemi == 'lh':
#                    offset = 0
#                    this_hemi = 0
#        
#                # get vertices on cortical surface inside label
#                idx = np.intersect1d(ll1.vertices, src_avg[this_hemi]['vertno'])
#        
#                # get vertices in source space inside label
#                fwd_idx = np.searchsorted(src_avg[this_hemi]['vertno'], idx)
#                thislabel_data[fwd_idx+offset,:]=datamat[llcnt,:].copy()
##                        lowsnr_thislabel_data[fwd_idx+offset,:]=lowsnr_datamat[llcnt,:].copy()
#            ## project from source space to sensor space
#            # multiply by leadfield
#            not_noise=np.unique(np.where(thislabel_data!=0)[0])
#            tobe_noise=np.asarray([thiszero for thiszero in range(thislabel_data.shape[0]) if thiszero not in not_noise])
#            tobe_allnoise=np.asarray([thiszero for thiszero in range(thislabel_data.shape[0])])
#
#            noise_std=this_std.copy()#np.sqrt(np.mean(datamat[0,blen:]**2))
##                    lowsnr_noise_std=lowsnr_this_std.copy()#np.sqrt(np.mean(nulldatamat[0,blen:]**2))
#            noise_mat=noise_std*np.random.randn(tobe_noise.shape[0],thislabel_data[:,blen:].shape[1])+np.mean(datamat[0,blen:])
##                    lowsnr_noise_mat=lowsnr_noise_std*np.random.randn(tobe_noise.shape[0],lowsnr_thislabel_data[:,blen:].shape[1])+np.mean(lowsnr_datamat[0,blen:])
#            
#            thislabel_data_ideal=thislabel_data.copy()
#            thislabel_data_ideal[tobe_noise,blen:]=0.0001*noise_mat
#            thislabel_data[tobe_noise,blen:]=noise_mat
#                
##                    lowsnr_thislabel_data[tobe_noise,blen:]=lowsnr_noise_mat
##                    allnoise_thislabel_data[:,blen:]=allnoise_mat
#            print ecnt
#            rmat=this_std*np.random.randn(thislabel_data.shape[0],thislabel_data.shape[1])
##                    lowsnr_rmat=lowsnr_this_std*np.random.randn(lowsnr_thislabel_data.shape[0],lowsnr_thislabel_data.shape[1])
#
#            thislabel_datan=thislabel_data+rmat
#            allnoise_std=np.mean(np.std(thislabel_datan[:,:blen],1))#baseline std
##                    allnoise_mn=np.mean(np.mean(thislabel_datan[:,:blen],1))#baseline std
#            allnoise_mat=allnoise_std*np.random.randn(tobe_allnoise.shape[0],thislabel_data.shape[1])#+np.mean(datamat[0,blen:])
#            allnoise_datan=allnoise_mat.copy()#allnoise_thislabel_data+rmat
#        
##                    lowsnr_thislabel_datan=lowsnr_thislabel_data+lowsnr_rmat
#            #morph created signals from avg to subj
#            # 1-high snr with noise
#            morphmat_avgtosub=mne.compute_morph_matrix(subject_from='fsaverage', subject_to=subjects[subject_no], vertices_from=vertices_avg, vertices_to=vertices_sub, subjects_dir=data_path)
#
#            stc_sim_pre=mne.SourceEstimate(thislabel_datan, vertices=vertices_avg, tmin=tmin, tstep=tstep, subject='fsaverage', verbose=None)#thislabel_datan
#            stc_sim_pre_morphed=stc_sim_pre.morph_precomputed(subject_to=subjects[subject_no],vertices_to=vertices_sub,morph_mat=morphmat_avgtosub, subject_from='fsaverage')
#            allnoise_stc_sim_pre=mne.SourceEstimate(allnoise_datan, vertices=vertices_avg, tmin=tmin, tstep=tstep, subject='fsaverage', verbose=None)#thislabel_datan
#            allnoise_stc_sim_pre_morphed=allnoise_stc_sim_pre.morph_precomputed(subject_to=subjects[subject_no],vertices_to=vertices_sub,morph_mat=morphmat_avgtosub, subject_from='fsaverage')
#
#            # 2- low snr with noise
##                    lowsnr_stc_sim_pre=mne.SourceEstimate(lowsnr_thislabel_datan, vertices=vertices_avg, tmin=tmin, tstep=tstep, subject='fsaverage', verbose=None)#thislabel_datan
##                    lowsnr_stc_sim_pre_morphed=lowsnr_stc_sim_pre.morph_precomputed(subject_to=subjects[subject_no],vertices_to=vertices_sub,morph_mat=morphmat_avgtosub, subject_from='fsaverage')
#            # 3- ideal 
##                    stc_sim_realidealpre=mne.SourceEstimate(thislabel_data_ideal, vertices=vertices_avg, tmin=tmin, tstep=tstep, subject=subjects[subject_no], verbose=None)#thislabel_datan
##                    stc_sim_realidealpre_morphed=stc_sim_realidealpre.morph_precomputed(subject_to=subjects[subject_no],vertices_to=vertices_sub,morph_mat=morphmat_avgtosub, subject_from='fsaverage')
#            stc_sim_ideal[ecnt]=mne.SourceEstimate(thislabel_datan, vertices=vertices_avg, tmin=tmin, tstep=tstep, subject='fsaverage', verbose=None)#thislabel_datan
#            stc_sim_realideal[ecnt]=mne.SourceEstimate(thislabel_data_ideal, vertices=vertices_avg, tmin=tmin, tstep=tstep, subject='fsaverage', verbose=None)#thislabel_datan
##                    lowsnr_stc_sim_ideal[ecnt]=mne.SourceEstimate(lowsnr_thislabel_datan, vertices=vertices_avg, tmin=tmin, tstep=tstep, subject='fsaverage', verbose=None)#thislabel_datan
#
#            thislabel_datan_morphed=stc_sim_pre_morphed.data.copy()
#            allnoise_thislabel_datan_morphed=allnoise_stc_sim_pre_morphed.data
#            evoked_data=np.matmul(leadfield,thislabel_datan_morphed)[np.newaxis,:,:]
##                    lowsnr_thislabel_datan_morphed=lowsnr_stc_sim_pre_morphed.data
##                    lowsnr_evoked_data=np.matmul(leadfield,lowsnr_thislabel_datan_morphed)[np.newaxis,:,:]
#            allnoise_evoked_data=np.matmul(leadfield,allnoise_thislabel_datan_morphed)[np.newaxis,:,:]
#            if ecnt==0:
#                epochs_data=evoked_data.copy()
##                        lowsnr_epochs_data=lowsnr_evoked_data.copy()
#                allnoise_epochs_data=allnoise_evoked_data.copy()
#                events=np.array([0,1,1])[np.newaxis,:]
#    
#            else:
#                epochs_data=np.append(epochs_data,evoked_data,axis=0)
#                allnoise_epochs_data=np.append(allnoise_epochs_data,allnoise_evoked_data,axis=0)
##                        lowsnr_epochs_data=np.append(lowsnr_epochs_data,lowsnr_evoked_data,axis=0)
#                this_event=np.array([ecnt+1,1,1])[np.newaxis,:]
#                events=np.append(events,this_event,axis=0)
#        # create evoked array
#    #    evoked_fwd = mne.EvokedArray(evoked_data, info=info, tmin=-0.769)
#        epochs_fwd = mne.EpochsArray(epochs_data, info=epochs1.info, events=events, event_id=1, tmin=tmin)
#        allnoise_epochs_fwd = mne.EpochsArray(allnoise_epochs_data, info=epochs1.info, events=events, event_id=1, tmin=tmin)
##                lowsnr_epochs_fwd = mne.EpochsArray(lowsnr_epochs_data, info=epochs1.info, events=events, event_id=1, tmin=tmin)
#
#        snr = thissnr#3.0#ediv.mean()
#        lambda2 = 1.0 / snr ** 2
##                lambda2lowsnr = 1.0 / lowsnr ** 2
#        method = 'MNE'
#        ## apply inverse to the evoked array
#        #stc_sim = apply_inverse(evoked_fwd, inverse_operator_eegmeg, lambda2,method=method, pick_ori="normal")
#        stc_sim = apply_inverse_epochs(epochs_fwd, inverse_operator_eegmeg, lambda2,method=method, pick_ori="normal")
#        allnoise_stc_sim = apply_inverse_epochs(allnoise_epochs_fwd, inverse_operator_eegmeg, lambda2,method=method, pick_ori="normal")
##                lowsnr_stc_sim = apply_inverse_epochs(lowsnr_epochs_fwd, inverse_operator_eegmeg, lambda2lowsnr,method=method, pick_ori="normal")
#        morphmat1=mne.compute_morph_matrix(subject_from=subjects[subject_no], subject_to='fsaverage', vertices_from=stc_sim[0].vertices, vertices_to=vertices_avg, subjects_dir=data_path)
#        stc_sim_morphed=[con1.copy().morph_precomputed(subject_to='fsaverage',vertices_to=vertices_avg,morph_mat=morphmat1, subject_from=subjects[subject_no]) for con1 in stc_sim]
#        allnoise_stc_sim_morphed=[con2.copy().morph_precomputed(subject_to='fsaverage',vertices_to=vertices_avg,morph_mat=morphmat1, subject_from=subjects[subject_no]) for con2 in allnoise_stc_sim]
#        del stc_sim
##                stc_sim_realideal_morphed=[con1.copy().morph_precomputed(subject_to='fsaverage',vertices_to=vertices_avg,morph_mat=morphmat1, subject_from=subjects[subject_no]) for con1 in stc_sim_realideal]
##                del stc_sim_realideal
##                lowsnr_stc_sim_morphed=[con1.copy().morph_precomputed(subject_to='fsaverage',vertices_to=vertices_avg,morph_mat=morphmat1, subject_from=subjects[subject_no]) for con1 in lowsnr_stc_sim]
##                del lowsnr_stc_sim
#        
#        print "extracting label time course"
#        label_ts = mne.extract_label_time_course(stc_sim_morphed, morphed_labelsall, src_avg, mode='mean_flip',return_generator=False)
#        allnoise_label_ts = mne.extract_label_time_course(allnoise_stc_sim_morphed, morphed_labelsall, src_avg, mode='mean_flip',return_generator=False)
##                lowsnr_label_ts = mne.extract_label_time_course(lowsnr_stc_sim_morphed, morphed_labelsall, src_avg, mode='mean_flip',return_generator=False)
##                lowsnr_label_ts_ideal = mne.extract_label_time_course(lowsnr_stc_sim_ideal, morphed_labelsall, src_sub, mode='mean_flip',return_generator=False)
##                label_ts_ideal = mne.extract_label_time_course(stc_sim_ideal, morphed_labelsall, src_sub, mode='mean_flip',return_generator=False)
#        label_ts_ideal = mne.extract_label_time_course(stc_sim_ideal, morphed_labelsall, src_avg, mode='mean_flip',return_generator=False)
#        label_ts_realideal = mne.extract_label_time_course(stc_sim_realideal, morphed_labelsall, src_avg, mode='mean_flip',return_generator=False)
##                lowsnr_label_ts_ideal = mne.extract_label_time_course(lowsnr_stc_sim_ideal, morphed_labelsall, src_avg, mode='mean_flip',return_generator=False)
#        
#            
#            ##################################
#            
#        #    cl=range(len(labels))
#        #    clw=range(len(labels))
#        #    cl_ideal=range(len(labels))
#        #    clw_ideal=range(len(labels))
#        #    cl_realideal=range(len(labels))
#        #    clw_realideal=range(len(labels))
#        #    for cntlts in range(n_epochs):
#        #        for cntcl in range(len(labels)):        
#        #            cl[cntcl]=simstc_morphed[cntlts].in_label(labels[cntcl]).data
#        #            clw[cntcl]=np.argmax(np.mean(np.abs(cl[cntcl][:,blen:]),axis=1))
#        #            thislabelsign=mne.label.label_sign_flip(labels[cntcl],src_avg)[clw[cntcl]]
#        #            label_ts[cntlts][cntcl,:]=thislabelsign*np.expand_dims(cl[cntcl][clw[cntcl],:].transpose(),0)
#        #            ##ideal
#        #            cl_ideal[cntcl]=simstc_morphed_ideal[cntlts].in_label(labels[cntcl]).data
#        #            clw_ideal[cntcl]=np.argmax(np.mean(np.abs(cl_ideal[cntcl][:,blen:]),axis=1))
#        #            thislabelsign2=mne.label.label_sign_flip(labels[cntcl],src_avg)[clw_ideal[cntcl]]
#        #            label_ts_ideal[cntlts][cntcl,:]=thislabelsign2*np.expand_dims(cl_ideal[cntcl][clw_ideal[cntcl],:].transpose(),0)
#        #            ##realideal
#        #            cl_realideal[cntcl]=simstc_morphed_realideal[cntlts].in_label(labels[cntcl]).data
#        #            clw_realideal[cntcl]=np.argmax(np.mean(np.abs(cl_realideal[cntcl][:,blen:]),axis=1))
#        #            thislabelsign2=mne.label.label_sign_flip(labels[cntcl],src_avg)[clw_realideal[cntcl]]
#        #            label_ts_realideal[cntlts][cntcl,:]=thislabelsign2*np.expand_dims(cl_realideal[cntcl][clw_realideal[cntcl],:].transpose(),0)
#
#        print "extracting label time course DONE"
#    
#        ##coherence
#        
#        fmin = 8.
#        fmax = 55.
#        sfreq = raw.info['sfreq']  # the sampling frequency
#        con_methods = ['coh']#['pli2_unbiased']
#        print "computing first coh"
#        
#        
#        nparcel=label_ts[0].shape[0]
#        conmine=np.zeros((nparcel,nparcel))
#        for thisroi1 in range(nparcel):
#            print thisroi1
#            for thisroi2 in range(thisroi1):#range(label_ts[0].shape[0]):
#                #if thisroi1>thisroi2:
#                label_ts_partial=list()
#                for thisepoch in range(len(label_ts)):
#                    x=label_ts[thisepoch].copy().transpose() #ntimesxnrois
#                    n,d=x.shape
#                    xx=x[:,[thisroi1,thisroi2]]
#                    zz=x[:,np.setdiff1d(np.arange(d),np.array([thisroi1,thisroi2]))]
#                    z1=np.hstack((np.ones((n,1)),zz))
#                    bb=np.linalg.lstsq(z1,xx)[0]
#                    resid=xx-np.matmul(z1,bb)
#                    dz=np.linalg.matrix_rank(x)-2
#                    tol=np.max([n,dz])*np.finfo(float).eps*np.sqrt(np.sum(np.abs(xx)**2,0))
#                    #tol1=np.tile(tol,(d,1))
#                    wzero=np.where(np.sqrt(np.sum(np.abs(resid)**2,0))<tol)[0]
#                    if len(wzero>0):
#                        resid[:,wzero]=0
#                    label_ts_partial.append(resid.copy().transpose())
#                thiscon, freqs, times, n_epochs, n_tapers = spectral_connectivity(
#                label_ts_partial, method=con_methods, mode='multitaper', sfreq=sfreq, fmin=fmin,fmax=fmax, 
#                faverage=True, mt_adaptive=True, n_jobs=1) 
#                conmine[thisroi1,thisroi2]=np.squeeze(thiscon[1,0,0].copy())
#                
#        con_noisemine=np.zeros((nparcel,nparcel))
#        for thisroi1 in range(nparcel):
#            print thisroi1
#            for thisroi2 in range(thisroi1):#range(label_ts[0].shape[0]):
#                #if thisroi1>thisroi2:
#                allnoise_label_ts_partial=list()
#                for thisepoch in range(len(allnoise_label_ts)):
#                    x=allnoise_label_ts[thisepoch].copy().transpose() #ntimesxnrois
#                    n,d=x.shape
#                    xx=x[:,[thisroi1,thisroi2]]
#                    zz=x[:,np.setdiff1d(np.arange(d),np.array([thisroi1,thisroi2]))]
#                    z1=np.hstack((np.ones((n,1)),zz))
#                    bb=np.linalg.lstsq(z1,xx)[0]
#                    resid=xx-np.matmul(z1,bb)
#                    dz=np.linalg.matrix_rank(x)-2
#                    tol=np.max([n,dz])*np.finfo(float).eps*np.sqrt(np.sum(np.abs(xx)**2,0))
#                    #tol1=np.tile(tol,(d,1))
#                    wzero=np.where(np.sqrt(np.sum(np.abs(resid)**2,0))<tol)[0]
#                    if len(wzero>0):
#                        resid[:,wzero]=0
#                    allnoise_label_ts_partial.append(resid.copy().transpose())
#                allnoise_thiscon, freqs, times, n_epochs, n_tapers = spectral_connectivity(
#                allnoise_label_ts_partial, method=con_methods, mode='multitaper', sfreq=sfreq, fmin=fmin,fmax=fmax, 
#                faverage=True, mt_adaptive=True, n_jobs=1) 
#                con_noisemine[thisroi1,thisroi2]=np.squeeze(allnoise_thiscon[1,0,0].copy())
#                
#        con_idealmine=np.zeros((nparcel,nparcel))
#        for thisroi1 in range(nparcel):
#            print thisroi1
#            for thisroi2 in range(thisroi1):#range(label_ts[0].shape[0]):
#                #if thisroi1>thisroi2:
#                ideal_label_ts_partial=list()
#                for thisepoch in range(len(label_ts_ideal)):
#                    x=label_ts_ideal[thisepoch].copy().transpose() #ntimesxnrois
#                    n,d=x.shape
#                    xx=x[:,[thisroi1,thisroi2]]
#                    zz=x[:,np.setdiff1d(np.arange(d),np.array([thisroi1,thisroi2]))]
#                    z1=np.hstack((np.ones((n,1)),zz))
#                    bb=np.linalg.lstsq(z1,xx)[0]
#                    resid=xx-np.matmul(z1,bb)
#                    dz=np.linalg.matrix_rank(x)-2
#                    tol=np.max([n,dz])*np.finfo(float).eps*np.sqrt(np.sum(np.abs(xx)**2,0))
#                    #tol1=np.tile(tol,(d,1))
#                    wzero=np.where(np.sqrt(np.sum(np.abs(resid)**2,0))<tol)[0]
#                    if len(wzero>0):
#                        resid[:,wzero]=0
#                    ideal_label_ts_partial.append(resid.copy().transpose())
#                ideal_thiscon, freqs, times, n_epochs, n_tapers = spectral_connectivity(
#                ideal_label_ts_partial, method=con_methods, mode='multitaper', sfreq=sfreq, fmin=fmin,fmax=fmax, 
#                faverage=True, mt_adaptive=True, n_jobs=1) 
#                con_idealmine[thisroi1,thisroi2]=np.squeeze(ideal_thiscon[1,0,0].copy())
#        
#            
#
#        con, freqs, times, n_epochs, n_tapers = spectral_connectivity(
#        label_ts, method=con_methods, mode='multitaper', sfreq=sfreq, fmin=fmin,fmax=fmax, 
#        faverage=True, mt_adaptive=True, n_jobs=1)                
#        con_res = np.squeeze(con.copy())#np.mean(con,axis=2)#dict()
#        print "computing second coh"
#        con_noise, freqs, times, n_epochs, n_tapers = spectral_connectivity(
#            allnoise_label_ts, method=con_methods, mode='multitaper', sfreq=sfreq, fmin=fmin,
#            fmax=fmax, faverage=True, mt_adaptive=True, n_jobs=1)
#        con_res_noise = np.squeeze(con_noise.copy())#np.mean(con_ideal,axis=2)#con_ideal[:,:,0]#dict()
#        print "computing 3rd coh"
#        con_ideal, freqs, times, n_epochs, n_tapers = spectral_connectivity(
#            label_ts_ideal, method=con_methods, mode='multitaper', sfreq=sfreq, fmin=fmin,
#            fmax=fmax, faverage=True, mt_adaptive=True, n_jobs=1)
#        con_res_ideal = np.squeeze(con_ideal.copy())
##                print "computing 4th coh"
##                lowsnr_con, freqs, times, n_epochs, n_tapers = spectral_connectivity(
##                lowsnr_label_ts, method=con_methods, mode='multitaper', sfreq=sfreq, fmin=fmin,fmax=fmax, 
##                faverage=True, mt_adaptive=True, n_jobs=1)                
##                lowsnr_con_res = np.squeeze(lowsnr_con.copy())
##                print "computing 5th coh"
##                lowsnr_con_ideal, freqs, times, n_epochs, n_tapers = spectral_connectivity(
##                lowsnr_label_ts_ideal, method=con_methods, mode='multitaper', sfreq=sfreq, fmin=fmin,fmax=fmax, 
##                faverage=True, mt_adaptive=True, n_jobs=1)                
##                lowsnr_con_res_ideal = np.squeeze(lowsnr_con_ideal.copy())#np.mean(con,axis=2)#dict()
#        
##                print "computing 1st imcoh"
##                con_methods2 = ['imcoh']#['pli2_unbiased']
##                con2,  freqs, times, n_epochs, n_tapers = spectral_connectivity(
##                    label_ts, method=con_methods2, mode='multitaper', sfreq=sfreq, fmin=fmin,
##                    fmax=fmax, faverage=True, mt_adaptive=True, n_jobs=1)
##            
##                con_res2 = np.squeeze(con2.copy())#np.mean(con2,axis=2)#con2[:,:,0]#dict()
##                print "computing second imcoh"
##                con_ideal2, freqs, times, n_epochs, n_tapers = spectral_connectivity(
##                    label_ts_ideal, method=con_methods2, mode='multitaper', sfreq=sfreq, fmin=fmin,
##                    fmax=fmax, faverage=True, mt_adaptive=True, n_jobs=1)
###                #     con is a 3D array, get the connectivity for the first (and only) freq. band
###                #     for each method
##                con_res_ideal2 = np.squeeze(con_ideal2.copy())#np.mean(con_ideal2,axis=2)#con_ideal2[:,:,0]#dict()
##                print "computing 3rd imcoh"
##                con_realideal2, freqs, times, n_epochs, n_tapers = spectral_connectivity(
##                    label_ts_ideal, method=con_methods2, mode='multitaper', sfreq=sfreq, fmin=fmin,
##                    fmax=fmax, faverage=True, mt_adaptive=True, n_jobs=1)
##                con_res_realideal2 = np.squeeze(con_realideal2.copy())    
##                print "computing 4th coh"
##                lowsnr_con2, freqs, times, n_epochs, n_tapers = spectral_connectivity(
##                lowsnr_label_ts, method=con_methods2, mode='multitaper', sfreq=sfreq, fmin=fmin,fmax=fmax, 
##                faverage=True, mt_adaptive=True, n_jobs=1)                
##                lowsnr_con_res2 = np.squeeze(lowsnr_con2.copy())
##                print "computing 5th coh"
##                lowsnr_con_ideal2, freqs, times, n_epochs, n_tapers = spectral_connectivity(
##                lowsnr_label_ts_ideal, method=con_methods2, mode='multitaper', sfreq=sfreq, fmin=fmin,fmax=fmax, 
##                faverage=True, mt_adaptive=True, n_jobs=1)                
##                lowsnr_con_res_ideal2 = np.squeeze(lowsnr_con_ideal2.copy())
##                #     con is a 3D array, get the connectivity for the first (and only) freq. band
##                #     for each method
#        
#        if which_parc==0:
#            parc0_connmat[:,:,0,cnt]=con_res.copy()#uv leakage
#            parc0_connmat[:,:,1,cnt]=con_res_noise.copy()#uv null
#            parc0_connmat[:,:,2,cnt]=con_res_ideal.copy()#uv ideal
#            parc0_connmat[:,:,3,cnt]=conmine.copy()#con_res_final=con_res[:,:,np.newaxis].copy()
#            parc0_connmat[:,:,4,cnt]=con_noisemine.copy()
#            parc0_connmat[:,:,5,cnt]=con_idealmine.copy()
#            parc0_connmat[:,:,6,cnt]=idealconn.copy()
#
#        if which_parc==1:
#            parc1_connmat[:,:,0,cnt]=con_res.copy()#con_res_final=con_res[:,:,np.newaxis].copy()
#            parc1_connmat[:,:,1,cnt]=con_res_noise.copy()
#            parc1_connmat[:,:,2,cnt]=con_res_ideal.copy()
#            parc1_connmat[:,:,3,cnt]=conmine.copy()#con_res_final=con_res[:,:,np.newaxis].copy()
#            parc1_connmat[:,:,4,cnt]=con_noisemine.copy()
#            parc1_connmat[:,:,5,cnt]=con_idealmine.copy()
#            parc1_connmat[:,:,6,cnt]=idealconn.copy()
#        if which_parc==2:
#            parc2_connmat[:,:,0,cnt]=con_res.copy()#con_res_final=con_res[:,:,np.newaxis].copy()
#            parc2_connmat[:,:,1,cnt]=con_res_noise.copy()
#            parc2_connmat[:,:,2,cnt]=con_res_ideal.copy()
#            parc2_connmat[:,:,3,cnt]=conmine.copy()#con_res_final=con_res[:,:,np.newaxis].copy()
#            parc2_connmat[:,:,4,cnt]=con_noisemine.copy()
#            parc2_connmat[:,:,5,cnt]=con_idealmine.copy()
#            parc2_connmat[:,:,6,cnt]=idealconn.copy()
#        if which_parc==3:
#            parc3_connmat[:,:,0,cnt]=con_res.copy()#con_res_final=con_res[:,:,np.newaxis].copy()
#            parc3_connmat[:,:,1,cnt]=con_res_noise.copy()
#            parc3_connmat[:,:,2,cnt]=con_res_ideal.copy()
#            parc3_connmat[:,:,3,cnt]=conmine.copy()#con_res_final=con_res[:,:,np.newaxis].copy()
#            parc3_connmat[:,:,4,cnt]=con_noisemine.copy()
#            parc3_connmat[:,:,5,cnt]=con_idealmine.copy()
#            parc3_connmat[:,:,6,cnt]=idealconn.copy()
#        if which_parc==4:
#            parc4_connmat[:,:,0,cnt]=con_res.copy()#con_res_final=con_res[:,:,np.newaxis].copy()
#            parc4_connmat[:,:,1,cnt]=con_res_noise.copy()
#            parc4_connmat[:,:,2,cnt]=con_res_ideal.copy()
#            parc4_connmat[:,:,3,cnt]=conmine.copy()#con_res_final=con_res[:,:,np.newaxis].copy()
#            parc4_connmat[:,:,4,cnt]=con_noisemine.copy()
#            parc4_connmat[:,:,5,cnt]=con_idealmine.copy()
#            parc4_connmat[:,:,6,cnt]=idealconn.copy()
#
##        this_pathto=out_path+outname[0]+'parc_connmat_sim'+str(perm_inds[0]+sim_offset)+'.mat'#'DKA68_orig_modified_coh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot            
##        scipy.io.savemat(this_pathto,{'parc0_connmat':parc0_connmat})
#this_pathto=out_path+outname[1]+'parc_connmat_sim'+str(perm_inds[0]+sim_offset)+'.mat'#'DKA68_orig_modified_coh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot            
#scipy.io.savemat(this_pathto,{'parc1_connmat':parc1_connmat})
##        this_pathto=out_path+outname[2]+'parc_connmat_sim'+str(perm_inds[0]+sim_offset)+'.mat'#'DKA68_orig_modified_coh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot            
##        scipy.io.savemat(this_pathto,{'parc2_connmat':parc2_connmat})
##        this_pathto=out_path+outname[3]+'parc_connmat_sim'+str(perm_inds[0]+sim_offset)+'.mat'#'DKA68_orig_modified_coh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot            
##        scipy.io.savemat(this_pathto,{'parc3_connmat':parc3_connmat})
##        this_pathto=out_path+outname[4]+'parc_connmat_sim'+str(perm_inds[0]+sim_offset)+'.mat'#'DKA68_orig_modified_coh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot            
##        scipy.io.savemat(this_pathto,{'parc4_connmat':parc4_connmat})
#import pickle
#this_pathto=out_path+'/activesources_grandmat_sim'+str(perm_inds[0]+sim_offset)+'.py'#'DKA68_orig_modified_coh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot            
#with open(this_pathto,'wb') as fp:
#    pickle.dump(actsrc_grandmat,fp)
#this_pathto=out_path+'/activesources_howmanyovs_sim'+str(perm_inds[0]+sim_offset)+'.py'#'DKA68_orig_modified_coh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot            
#with open(this_pathto,'wb') as fp:
#    pickle.dump(actsrc_ovs,fp)        
#this_pathto=out_path+'/activesourcesfinal_grandmat_sim'+str(perm_inds[0]+sim_offset)+'.py'#'DKA68_orig_modified_coh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot            
#with open(this_pathto,'wb') as fp:
#    pickle.dump(actsrcfinal_grandmat,fp)
#this_pathto=out_path+'/missedcons_grandmat_sim'+str(perm_inds[0]+sim_offset)+'.py'#'DKA68_orig_modified_coh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot            
#with open(this_pathto,'wb') as fp:
#    pickle.dump(missedcons_grandmat,fp)
#missedcons_gmat=np.asarray(missedcons_grandmat)   
#this_pathto=out_path+'/missedcons_grandmat_sim'+str(perm_inds[0]+sim_offset)+'.mat'#'DKA68_orig_modified_coh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot            
#scipy.io.savemat(this_pathto,{'missedcons_gmat':missedcons_gmat})
#this_pathto=out_path+'/misscatch_grandmat_sim'+str(perm_inds[0]+sim_offset)+'.py'#'DKA68_orig_modified_coh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot            
#with open(this_pathto,'wb') as fp:
#    pickle.dump(misscatch_grandmat,fp)
#
#    
#
#  
##with open ('outfile', 'rb') as fp:
##    itemlist = pickle.load(fp)
## First, we reorder the labels based on their location in the left hemi
##tfin=time.time()-t0
##label_names = [label.name for label in labels]
##
##lh_labels = [name for name in label_names if name.endswith('lh')]
##
### Get the y-location of the label
##label_ypos = list()
##for name in lh_labels:
##    idx = label_names.index(name)
##    ypos = np.mean(labels[idx].pos[:, 1])
##    label_ypos.append(ypos)
##
### Reorder the labels based on their location
##lh_labels = [label for (yp, label) in sorted(zip(label_ypos, lh_labels))]
##
### For the right hemi
##rh_labels = [label[:-2] + 'rh' for label in lh_labels]
##
### Save the plot order and create a circular layout
##node_order = list()
##node_order.extend(lh_labels[::-1])  # reverse the order
##node_order.extend(rh_labels)
##
##node_angles = circular_layout(label_names, node_order, start_pos=90,
##                              group_boundaries=[0, len(label_names) / 2])
##
### Plot the graph using node colors from the FreeSurfer parcellation. We only
### show the 300 strongest connections.
##plt.close("all")
##con_res_mean=np.squeeze(np.mean(con_res_final,axis=2))
##plot_connectivity_circle(np.abs(con_res_mean), range(len(label_names)),n_lines=None,
##                         node_angles=node_angles, node_colors=label_colors,#vmin=0,vmax=1,
##                         title='All-to-All Connectivity '
##                               'Coh')
###con_mod_mean=np.mean(con_mod_final,axis=2)                             
###plot_connectivity_circle(np.abs(con_mod_mean), range(len(label_names)),n_lines=300,
###                         node_angles=node_angles, node_colors=label_colors,#vmin=0,vmax=1e-19,
###                         title='All-to-All Connectivity '
###                               'Coh')
##con_res_mean_ideal=np.squeeze(np.mean(con_res_final_ideal,axis=2))
##plot_connectivity_circle(np.abs(con_res_mean_ideal), label_names, n_lines=None,#vmin=0,#vmax=1,
##                         node_angles=node_angles, node_colors=label_colors,
##                         title='All-to-All Connectivity no leakage Coh')
#####
##con_res_mean_realideal=np.squeeze(np.mean(con_res_final_realideal,axis=2))
##plot_connectivity_circle(np.abs(con_res_mean_realideal), label_names, n_lines=None,#vmin=0,#vmax=1,
##                         node_angles=node_angles, node_colors=label_colors,
##                         title='All-to-All Connectivity '
##                               'Ideal Coh')                               
#####
##con_res_mean2=np.squeeze(np.mean(con_res_final2,axis=2))
##plot_connectivity_circle(np.abs(con_res_mean2), label_names, n_lines=None,
##                         node_angles=node_angles, node_colors=label_colors,#vmin=0,vmax=1,
##                         title='All-to-All Connectivity '
##                               'ImCoh')
###con2_mod_mean=np.mean(con2_mod_final,axis=2)                             
###plot_connectivity_circle(np.abs(con2_mod_mean), range(len(label_names)),n_lines=300,
###                         node_angles=node_angles, node_colors=label_colors,#vmin=0,vmax=1,
###                         title='All-to-All Connectivity '
###                               'Coh')
##con_res_mean_ideal2=np.squeeze(np.mean(con_res_final_ideal2,axis=2))
##plot_connectivity_circle(np.abs(con_res_mean_ideal2), label_names, n_lines=None,#vmin=0,vmax=1,
##                         node_angles=node_angles, node_colors=label_colors,
##                         title='All-to-All Connectivity '
##                               'no leakage ImCoh') 
##                               
#####
##con_res_mean_realideal2=np.squeeze(np.mean(con_res_final_realideal2,axis=2))
##plot_connectivity_circle(np.abs(con_res_mean_realideal2), label_names, n_lines=None,#vmin=0,vmax=1,
##                         node_angles=node_angles, node_colors=label_colors,
##                         title='All-to-All Connectivity '
##                               'Ideal ImCoh')                                
###plt.savefig('circle.png', facecolor='black')
##outname=['aparc','aparc_mod','aparc2009','aparc2009_mod','RG']
###out_path='/imaging/rf02/parcellation/extended_sim/'
###path_to=label_path+outname[which_parc]+'_coh_test.mat'#'DKA68_orig_modified_coh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot
###scipy.io.savemat(path_to,{'con_res_final':con_res_final})
####path_to=label_path+outname[which_parc]+'_coh_mod_test.mat'#'DKA68_orig_modified_coh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot
####scipy.io.savemat(path_to,{'con_mod_final':con_mod_final})
###
###path_to=label_path+outname[which_parc]+'_ideal_coh_test.mat'#'DKA68_ideal_modified_coh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot
###scipy.io.savemat(path_to,{'con_res_final_ideal':con_res_final_ideal})
###
###path_to=label_path+outname[which_parc]+'_imcoh_test.mat'#'DKA68_orig_modified_imcoh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot
###scipy.io.savemat(path_to,{'con_res_final2':con_res_final2})
####path_to=label_path+outname[which_parc]+'_imcoh_mod_test.mat'#'DKA68_orig_modified_imcoh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot
####scipy.io.savemat(path_to,{'con2_mod_final':con2_mod_final})
###
###path_to=label_path+outname[which_parc]+'_ideal_imcoh_test.mat'#'DKA68_ideal_modified_imcoh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot
###scipy.io.savemat(path_to,{'con_res_final_ideal2':con_res_final_ideal2})
##plt.show()
## 
##def mean2(x):
##    y=np.sum(x)/np.size(x)
##    return y
##def corr2(a,b):
##    a=a-mean2(a)
##    b=b-mean2(b)
##    r=(a*b).sum()/np.sqrt((a*a).sum()*(b*b).sum())
##    return r       
##cc2=(con_res_mean_realideal+con_res_mean_realideal.T)/2
###cc2=cc2+np.eye(cc2.shape[0])
##cc1=(con_res_mean+con_res_mean.T)/2 
###cc1=cc1+np.eye(cc1.shape[0])  
##R=corr2(cc1,cc2)
#  #                xx=x[:,[ii,jj]]
##                zz=x[:,:-2]
##                    xx.shape
##(100, 2)
##>>> zz.shape
##(100, 18)
##>>> z1=np.hstack((np.ones((d,1)),zz))
##>>> z1.shape
##(100, 19)
##>>> bb=np.linalg.lstsq(z1,xx)[0]