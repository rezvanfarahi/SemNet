"""
=========================================================
TF representation of SemLoc for frequency bands in source labels
=========================================================

"""
# Authors: Alexandre Gramfort <gramfort@nmr.mgh.harvard.edu>
#
# License: BSD (3-clause)


print(__doc__)

import sys
sys.path.append('/imaging/local/software/python_packages/nibabel/1.3.0')
sys.path.append('/imaging/local/software/python_packages/pysurfer/v0.3.1')
sys.path.append('/imaging/local/software/python_packages/scikit-learn/v0.14.1')
sys.path.insert(1,'/imaging/rf02/TypLexMEG/mne_python0.9')

# for qsub
# (add ! here if needed) 
sys.path.append('/imaging/local/software/anaconda/latest/x86_64/bin/python')
sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
sys.path.insert(1,'/imaging/local/software/anaconda/2.4.1/3/lib/python3.5/site-packages/scipy')
###

sys.path.insert(1,'/imaging/local/software/python_packages/scikit-learn/v0.13.1/lib.linux-x86_64-2.7/sklearn/externals')
import joblib
###


import os
import numpy as np
import mne
import scipy
from scipy.stats.mstats import zscore
import scipy.io as scio
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
import copy

###############################################################################
out_path = '/home/rf02/rezvan/test1/step_by_step/dcm/' # root directory for your MEG data
subjects_dir = '/imaging/rf02/TypLexMEG/'    # where your MRI subdirectories are
# where event-files aremean
data_path = '/imaging/rf02/TypLexMEG/'    # where event files are

label_path = '/home/rf02/rezvan/TypLexMEG/createdlabels_SL/'

inv_path = '/imaging/rf02/TypLexMEG/'

# get indices for subjects to be processed from command line input
# 
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
for ss in subject_inds:
 ll.append(list_all[ss])

print "ll:"
print ll
labelss = mne.read_labels_from_annot(subject='fsaverage', parc='rez_aparc_25oct16', subjects_dir=data_path)
#labelss.remove(labelss[-3])
#labelss.remove(labelss[-2])
for lbl in labelss:
    lbl.values.fill(1.0)
label_names=[label.name for label in labelss]
label_index=np.array([label_names.index('ATL_latmed-lh'),label_names.index('posterior_inferiorfrontal-lh'),label_names.index('middle_temporal-lh'),label_names.index('AG_new-lh'),label_names.index('vWFA2-lh')])#label_names.index('ATL_latmed-lh'),label_names.index('supramarginal-lh'),

labellist=[labelss[ii] for ii in label_index]
#labellist = ['lh.ATL_latmed.label']#['atlleft-lh.label']#['lh.ATL_latmed.label']
#labellist2 = ['AG_new-lh.label','vWFA2-lh.label']#'AGSMG_meta2-lh.label'['lh.ATL_DCM-lh.label','lh.WFA_DCM-lh.label']
#lfname=data_path+'fsaverage/label/'+labellist2[0]#label_path+labellist[0]#data_path+'fsaverage/label/'+labellist[0]
#label0=mne.read_label(lfname,subject='fsaverage')
#labellist.append(label0)
#lfname=data_path+'fsaverage/label/'+labellist2[1]#label_path+labellist[0]#data_path+'fsaverage/label/'+labellist[0]
#label1=mne.read_label(lfname,subject='fsaverage')
#labellist.append(label1)
#print labellist
#labellist = ['lh.WFA_DCM-lh.label','lh.ATL_DCM-lh.label']

# Labels/ROIs
#labellist = ['atlleft-lh']#, 'atlleftt-rh'
tmin, tmax = -0.5, 0.7
# frequency bands with variable number of cycles for wavelets
stim_delay = 0.034 # delay in s

tmin1=-500
tstep1=700
stc_allc=range(17)
vertices_to = [np.arange(10242), np.arange(10242)]
stc_alla=range(17)
vertices_to = [np.arange(10242), np.arange(10242)]
cl=range(5)
al=range(5)
clw=np.arange(5)
alw=np.arange(5)
mlw=np.arange(5)
cncrt_mat=np.zeros((5,1201))
labels_tcc=np.zeros((17,1201,5))
labels_tca=np.zeros((17,1201,5))

cons=np.zeros((17,5,5))
abs_mat=np.zeros((5,1201))
conabs_mat=np.zeros((5,1201,2))
sfrq=1.
n_subjects=17
n_times=int(1201/sfrq)
n_levels_facts=4
X=np.zeros((n_subjects,5,n_times,2))
sonset=int(abs(tmin*1000)/sfrq)
import matplotlib.pyplot as plt

for ii, meg in enumerate(ll):
    conabs_mat=np.zeros((5,1201,2))
    subject_no=subject_inds[ii]
#    conind=0
    for event_no,event_name in enumerate(['Concrete','Abstract']):
#        conind=conind+1
#        print conind
        snr = 3.0
        print snr
        lambda2 = 1.0 / snr ** 2
        fname_fwd = data_path + meg + 'ico5_forward_5-3L-EMEG-fwd.fif'
        forward = mne.read_forward_solution(fname_fwd,surf_ori=True, force_fixed=True)
        print "load source estimates"
        method = "MNE"  # use dSPM method (could also be MNE or sLORETA)
        fname_inv = data_path + meg + 'InvOp_ico5diagreg_fft_1_48_clean_ica_EMEG-inv.fif'#'InvOp_ico5diagreg_fft_1_48_clean_ica_EMEG-inv.fif'#
        # Load data
        inverse_operator = read_inverse_operator(fname_inv)
        srcin='/imaging/rf02/TypLexMEG/fsaverage/bem/fsaverage_dist-ico-5-src.fif'
        src = mne.read_source_spaces(srcin)
        fname = data_path + meg + 'firstMorphed_ico_signed_SemLoc_ica_'+event_name+'_Source_Evoked_m500_700' 
        this_stc = mne.read_source_estimate(fname)  
        tmin1=-500
        tstep1=1
        thisstc_data=np.ones((this_stc.data.shape[0],1))
        vertices_to = [np.arange(10242), np.arange(10242)]
        thisstc_tolabel = mne.SourceEstimate(thisstc_data, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
        thisstc_label=mne.stc_to_label(thisstc_tolabel, src)

        ntimes=this_stc.data.shape[1]
        print "select 100 most sensitive labels"
        
        method = 'MNE'  # can be 'MNE' or 'sLORETA'
        mode = 'svd'
        n_svd_comp = 1
        labellist_temp=copy.deepcopy(labellist)
        labels = [labellist_temp[jj].morph(subject_from='fsaverage', subject_to=subjects[subject_no], subjects_dir=data_path) for jj in range(len(labellist))]
        stc_ctf_mne=psf_ctf_28Oct15.cross_talk_function(inverse_operator, forward, labels,method=method, lambda2=lambda2, signed=False, mode=mode, n_svd_comp=1, verbose=None)
        print labellist
        
        subject_from = subjects[subject_no]
        subject_to = 'fsaverage'
        vertices_to = [np.arange(10242), np.arange(10242)]
        morphmat1=mne.compute_morph_matrix(subject_from, subject_to, stc_ctf_mne.vertices, vertices_to, subjects_dir=data_path)
        stc_to=stc_ctf_mne.morph_precomputed(subject_to,vertices_to,morphmat1, subject_from)
        Matx=stc_to.data[:,:-1]# Nv x Nlabels+1
        srcin='/imaging/rf02/TypLexMEG/fsaverage/bem/fsaverage_dist-ico-5-src.fif'
        src_avg = mne.read_source_spaces(srcin)
        
        for this_labelc, this_label in enumerate(labellist): 
            label_origdata=this_stc.in_label(this_label).data
            label_verts=this_stc.in_label(this_label).lh_vertno
#            np.intersect1d(label_verts,vertices_to[0])
            vert25=np.sort(Matx[label_verts,this_labelc])[-25]
            winner_verts=np.where(Matx[:,this_labelc]>=np.max(Matx[:,this_labelc])/2.)[0]#label_verts[np.where(Matx[label_verts,this_labelc]>=np.max(Matx[label_verts,this_labelc])/2.)[0]]#label_verts[np.where(Matx[label_verts,this_labelc]>=vert25)[0]]#
            print winner_verts
            thisscale=np.hstack((mne.label.label_sign_flip(thisstc_label[0],src_avg),mne.label.label_sign_flip(thisstc_label[1],src_avg)))[winner_verts]#mne.label.label_sign_flip(this_label,src_avg)[winner_verts]#1#Matx[winner_verts,this_labelc]/np.max(Matx[winner_verts,this_labelc])
            label_predata=this_stc.data[winner_verts,:]
#            label_data1=np.tile(thisscale,[ntimes,1]).transpose()*label_predata
            label_data1=thisscale[:,np.newaxis]*label_predata
            active_ones=np.where(np.mean(np.abs(label_data1[:,sonset:]),1)>=np.mean(np.mean(np.abs(label_data1[:,sonset:]),1)))[0]
            label_data=label_data1[active_ones,:]
            ini_tc=np.mean(label_data[:,:], axis=0)
            plt.figure()
            plt.plot(ini_tc)
            conabs_mat[this_labelc,:,event_no]=ini_tc#(ini_tc-np.min(ini_tc))/(np.max(ini_tc)-np.min(ini_tc))zscore(ini_tc)#
            thisstc_testdata=np.zeros((this_stc.data.shape[0],1))
            thisstc_testdata[winner_verts[active_ones],0]=1
            vertices_test = [np.arange(10242), np.arange(10242)]
            thisstc_tolabel_test = mne.SourceEstimate(thisstc_testdata, vertices=vertices_test,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
            thisstc_labeltest=mne.stc_to_label(thisstc_tolabel_test, src_avg,smooth=False)[0]
            testlabel_path=data_path+'dcm/testDCMoldctflabs/'+'subject'+str(ii)+'_'+event_name+'_'+this_label.name
            thisstc_labeltest.save(testlabel_path)
            
#        wnw_mat=np.concatenate((np.expand_dims(cncrt_mat,2),np.expand_dims(abs_mat,2)),2)
#        con=mne.connectivity.spectral_connectivity(data=cncrt_mat[np.newaxis,:,550:750],method='coh',sfreq=1000.,mode='multitaper',fmin=1.,fmax=45.,mt_adaptive=True,faverage=True)
#        conc=np.corrcoef(cncrt_mat[:,550:750])#con[0]+con[0].transpose(1,0,2)
#        con2=mne.connectivity.spectral_connectivity(data=abs_mat[np.newaxis,:,550:750],method='coh',sfreq=1000.,mode='multitaper',fmin=1.,fmax=45.,mt_adaptive=True,faverage=True)
#        cona=np.corrcoef(abs_mat[:,550:750])#con2[0]+con2[0].transpose(1,0,2)
#        cons[ii,:,:]=np.squeeze(conc-cona)
#    X[ii,:,:,:]=conabs_mat
#    out_file=out_path + meg[:10] + '_SemLoc_Evoked_5ROIs_meanCTF_forDCM_new.mat'
#    scio.savemat(out_file,{'conabs_mat':conabs_mat})
#xms=np.mean(X,0)
#xms=np.mean(xms,2)
#import matplotlib.pyplot as plt
#plt.figure(1),plt.plot(np.linspace(-500,700,1201),xms[0,:])
#plt.figure(2),plt.plot(np.linspace(-500,700,1201),xms[1,:])
#plt.figure(3),plt.plot(np.linspace(-500,700,1201),xms[2,:])
#plt.figure(4),plt.plot(np.linspace(-500,700,1201),xms[3,:])
#plt.figure(5),plt.plot(np.linspace(-500,700,1201),xms[4,:])
#            conind=conind+1
#cl_ideal[cntcl]=simstc_morphed_ideal[cntlts].in_label(labels[cntcl]).data
#clw_ideal[cntcl]=np.argmax(np.mean(np.abs(cl_ideal[cntcl][:,125:]),axis=1))
#thislabelsign2=mne.label.label_sign_flip(labels[cntcl],src_avg)[clw_ideal[cntcl]]
#label_ts_ideal[cntlts][cntcl,:]=thislabelsign2*np.expand_dims(cl_ideal[cntcl][clw_ideal[cntcl],:].transpose(),0)

#    print ii
#    srcin='/imaging/rf02/TypLexMEG/fsaverage/bem/fsaverage-ico-5-src.fif'
#    src_avg = mne.read_source_spaces(srcin)
#    fname = data_path + meg + 'firstMorphed_ico_signed_SemLoc_ica_Concrete_Source_Evoked_m500_700' 
#    stc_allc = mne.read_source_estimate(fname)
#    thislabel_tcc=mne.extract_label_time_course(stc_allc, thislabel,src, mode='mean_flip') 
#    labels_tcc[ii,:,thislabelc]=thislabel_tcc.copy()
#    
##        cl1=stc_allc.in_label(thislabel).data
##    cl1w=np.argmax(np.mean(np.abs(cl1[:,550:950]),axis=1))
##    print cl1w
#    
#    
#    fname = data_path + meg + 'firstMorphed_ico_signed_SemLoc_ica_Abstract_Source_Evoked_m500_700'  
#    stc_alla = mne.read_source_estimate(fname)
#    thislabel_tca=mne.extract_label_time_course(stc_alla, thislabel,src, mode='mean_flip')
#    labels_tca[ii,:,thislabelc]=thislabel_tca.copy()
#
#labelsa_forsign=np.squeeze(np.mean(labels_tca,axis=0)) #size time x ROI
#labelsc_forsign=np.squeeze(np.mean(labels_tcc,axis=0)) #size time x ROI
#    
#for thislabelc,thislabel in enumerate(labellist):
#    print thislabel
#    for ii, meg in enumerate(ll):
#        print ii
#        srcin='/imaging/rf02/TypLexMEG/fsaverage/bem/fsaverage-ico-5-src.fif'
#        src = mne.read_source_spaces(srcin)
#        fname = data_path + meg + 'firstMorphed_ico_signed_SemLoc_ica_Concrete_Source_Evoked_m500_700' 
#        stc_allc = mne.read_source_estimate(fname)
#        thislabel_tcc=mne.extract_label_time_course(stc_allc, thislabel,src, mode='pca_flip') 
#        #        thislabel_signtc=mne.extract_label_time_course(stc_allc, thislabel,src, mode='mean_flip')
#        Pc=np.polyfit(np.linspace(50,200,150),thislabel_tcc[0,550:700],1)[0]
#        #        Pct=np.polyfit(np.linspace(50,200,150),thislabel_signtc[0,550:700],1)[0]
#        if np.corrcoef(np.squeeze(thislabel_tcc),labelsc_forsign[:,thislabelc])[0,1]>=0:#np.corrcoef(thislabel_tcc,thislabel_signtc)[0,1]>=0:#np.sign(Pc)==np.sign(Pct):
#            cncrt_mat[thislabelc,:]=thislabel_tcc.copy()
#            labels_tcc[ii,:,thislabelc]=thislabel_tcc.copy()
#        else:
#            cncrt_mat[thislabelc,:]=-thislabel_tcc.copy()
#            labels_tcc[ii,:,thislabelc]=-thislabel_tcc.copy()
#
#        #        cl1=stc_allc.in_label(thislabel).data
#        #    cl1w=np.argmax(np.mean(np.abs(cl1[:,550:950]),axis=1))
#        #    print cl1w
#        
#        
#        fname = data_path + meg + 'firstMorphed_ico_signed_SemLoc_ica_Abstract_Source_Evoked_m500_700'  
#        stc_alla = mne.read_source_estimate(fname)
#        thislabel_tca=mne.extract_label_time_course(stc_alla, thislabel,src, mode='pca_flip')
#        #        thislabel_signta=mne.extract_label_time_course(stc_alla, thislabel,src, mode='mean_flip')
#        Pa=np.polyfit(np.linspace(50,200,150),thislabel_tca[0,550:700],2)[0]
#        #        Pat=np.polyfit(np.linspace(50,200,150),thislabel_signta[0,550:700],2)[0]
#        if np.corrcoef(np.squeeze(thislabel_tca),labelsa_forsign[:,thislabelc])[0,1]>=0:#np.corrcoef(thislabel_tca,thislabel_signta)[0,1]>=0:#np.sign(Pa)==np.sign(Pat):
#            labels_tca[ii,:,thislabelc]=thislabel_tca.copy()
#            abs_mat[thislabelc,:]=thislabel_tca.copy()
#        else:
#            labels_tca[ii,:,thislabelc]=-thislabel_tca.copy()
#            abs_mat[thislabelc,:]=-thislabel_tca.copy()
        
#        wnw_mat=np.concatenate((np.expand_dims(cncrt_mat,2),np.expand_dims(abs_mat,2)),2)
#        con=mne.connectivity.spectral_connectivity(data=cncrt_mat[np.newaxis,:,550:750],method='coh',sfreq=1000.,mode='multitaper',fmin=1.,fmax=45.,mt_adaptive=True,faverage=True)
#        conc=np.corrcoef(cncrt_mat[:,550:750])#con[0]+con[0].transpose(1,0,2)
#        con2=mne.connectivity.spectral_connectivity(data=abs_mat[np.newaxis,:,550:750],method='coh',sfreq=1000.,mode='multitaper',fmin=1.,fmax=45.,mt_adaptive=True,faverage=True)
#        cona=np.corrcoef(abs_mat[:,550:750])#con2[0]+con2[0].transpose(1,0,2)
#        cons[ii,:,:]=np.squeeze(conc-cona)
#        
#        out_file=out_path + meg[:10] + '_SemLoc_Evoked_5ROIs_PCA_forDCM.mat'
#        scio.savemat(out_file,{'wnw_mat':wnw_mat})

#consr=np.reshape(cons,(17,25))
#import scipy.stats as scist
#MTc= mne.stats.ttest_1samp_no_p(consr, sigma=1e-3, method='relative')#Tc,pTc,tTc= mne.stats.permutation_t_test(consr)
#pTc=scist.t.sf(np.abs(MTc),16)*2
#pfinal=np.reshape(pTc,(5,5))
#print pfinal
#import matplotlib.pyplot as plt
#plt.figure(3)
#for ii in range(17):
#    plt.subplot(4,5,ii+1) 
#    plt.plot(np.linspace(-500,700,1201),labels_tcc[ii,:,4].T)
#    
#plt.figure(4)
#for ii in range(17):
#    plt.subplot(4,5,ii+1) 
#    plt.plot(np.linspace(-500,700,1201),labels_tca[ii,:,4].T)