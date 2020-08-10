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
#sys.path.insert(1,'/imaging/rf02/mne_python_11')
#sys.path.insert(1,'/imaging/local/software/mne_python/la_v0.9full')
sys.path.insert(1,'/imaging/local/software/mne_python/latest_v0.13')


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

import matplotlib.pyplot as plt
import copy
import pickle
print mne.__path__
import time
import math
def make_parcellation_mine(src,labels,colors,subject,fname):
    annotl=np.zeros((src[0]['np'],1))
    annotr=np.zeros((src[1]['np'],1))
    CTAB_L=np.zeros((len(labels)/2,5))
    CTAB_R=np.zeros((len(labels)/2,5))
    annot_id_coding = np.array((1, 2 ** 8, 2 ** 16))
    annot_id=np.sum(colors[:,:3] * annot_id_coding, axis=1)
    for lcntn in range(len(labels)):
        print lcntn
        if labels[lcntn].name.endswith('lh'):
            annotl[labels[lcntn].vertices]=annot_id[lcntn/2]
        else:
            annotr[labels[lcntn].vertices]=annot_id[int(lcntn/2)]
    CTAB_L[:,:4]=colors; CTAB_L[:,4]=annot_id
    CTAB_R[:,:4]=colors; CTAB_R[:,4]=annot_id
    label_names=[labels[lh].name[:-3] for lh in range(len(labels))]
    NAMES_L=label_names[::2]
    NAMES_R=label_names[1::2]
    ## write morphed labels
    fnamel=data_path + subject+'/label/lh.'+subject+'_'+fname+'.annot'
    mnewa(fnamel, np.squeeze(annotl), CTAB_L, NAMES_L)
    
    fnamer=data_path + subject+'/label/rh.'+subject+'_'+fname+'.annot'
    mnewa(fnamer, np.squeeze(annotr), CTAB_R, NAMES_R)

data_path = '/imaging/rf02/TypLexMEG/' # root directory for your MEG data
os.chdir(data_path)
event_path = '/group/erp/data/olaf.hauk/MEG/TypLexM/data/MF22_new/'
#label_path = '/group/erp/data/olaf.hauk/MEG/TypLexM/labels/createdlabels_SL/'
inv_path = '/imaging/rf02/TypLexMEG/'
subjects_dir = data_path 
print sys.argv
perm_inds = []
for ss in sys.argv[1:]:
   perm_inds.append( int( ss ) )

subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
print "subject_inds:"
print subject_inds
print "No rejection"
t0=time.time()
n_sims=70
sim_offset=0
sim_all=list(np.asarray(range(n_sims))+sim_offset)
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
for ss in perm_inds:
 ll.append(sim_all[ss])


print "ll:"
print ll

label_path='/imaging/rf02/TypLexMEG/parcellation/labels_and_CTFs/'
results_path='/imaging/rf02/parcellation/leakage_pattern/'

##grow ROIs
srcin='/imaging/rf02/TypLexMEG/fsaverage/bem/fsaverage_dist-ico-5-src.fif'
#src_avg2=src_avg.copy()
src_avg = mne.read_source_spaces(srcin)
subject_avg='fsaverage'
#nseeds=20
ext=25
connperc=0.25#0.25


#labellist_path = [label_path+'rostfront_smallseed-lh.label',label_path+'midtemp_smallseed-lh.label',label_path+'latocci_smallseed-lh.label']#[ 'lh.posterior_fusiform','rh.posterior_fusiform', 'lh.lateraloccipital','rh.lateraloccipital','lh.precentral','rh.precentral', 'lh.ATL_latmed','rh.ATL_latmed']#labellist = [ 'leftATL-lh','rightATL-rh']#, 'medialtempright-rh','medialtempleft-lh']'lh.postcentral','rh.postcentral',
#labellist_path = [label_path+'IFG_smallseed-lh.label',label_path+'IPL_smallseed-lh.label']#[ 'lh.posterior_fusiform','rh.posterior_fusiform', 'lh.lateraloccipital','rh.lateraloccipital','lh.precentral','rh.precentral', 'lh.ATL_latmed','rh.ATL_latmed']#labellist = [ 'leftATL-lh','rightATL-rh']#, 'medialtempright-rh','medialtempleft-lh']'lh.postcentral','rh.postcentral',
#labellist_path = [label_path+'LO_smallseed-lh.label',label_path+'RMF_smallseed-lh.label',label_path+'STS_smallseed-lh.label']#[ 'lh.posterior_fusiform','rh.posterior_fusiform', 'lh.lateraloccipital','rh.lateraloccipital','lh.precentral','rh.precentral', 'lh.ATL_latmed','rh.ATL_latmed']#labellist = [ 'leftATL-lh','rightATL-rh']#, 'medialtempright-rh','medialtempleft-lh']'lh.postcentral','rh.postcentral',
#labellist=[mne.read_label(label) for label in labellist_path]
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
## preparing unknown
labels_aparc = mne.read_labels_from_annot('fsaverage', parc=parc_name[0], subjects_dir=subjects_dir)
aparc_all=np.sum(labels_aparc)
aparc_l=aparc_all.lh
idx = np.intersect1d(aparc_l.vertices, src_avg[0]['vertno'])
fwd_idxl = np.searchsorted(src_avg[0]['vertno'], idx)
aparc_r=aparc_all.rh
offset = src_avg[0]['vertno'].shape[0]
idx = np.intersect1d(aparc_r.vertices, src_avg[1]['vertno'])
fwd_idxr = np.searchsorted(src_avg[1]['vertno'], idx)
#fwd_idxr=fwd_idxrp+offset
# get vertices on cortical surface inside label




parc_len=np.zeros((len(parc_name),))
for which_parc in range(len(parc_name)):
    labels_len = mne.read_labels_from_annot('fsaverage', parc=parc_name[which_parc], subjects_dir=subjects_dir)
    if which_parc==0:
        labels_len.remove(labels_len[-1])
    if which_parc==2:
        labels_len.remove(labels_len[-1])
        labels_len.remove(labels_len[-1])
    parc_len[which_parc]=len(labels_len)
    

            
            
            
            
#labels_pre.remove(labels_pre[-1])
#    labels_pre.remove(labels_pre[-1])

for simcnt in ll:
    for nseeds in np.array([15]): #,10,15,20
        print nseeds
#        thissnr=3
        out_path='/imaging/rf02/parcellation/extended_simulations/sims'+str(n_sims)+'_offset'+str(sim_offset)+'_seeds'+str(nseeds)+'_conns'+str(connperc)+'_bothsnrs_ext'+str(ext)
        outname=['/aparc/','/aparc_mod/','/aparc2009/','/aparc2009_mod/','/RG/']
        for thisname in outname:
            if not os.path.exists(out_path+thisname):
                os.makedirs(out_path+thisname)
        #mne.write_labels_to_annot(morphed_labelsall,subject=subjects[subject_no],parc='test_parc',subjects_dir=data_path,overwrite=True)
        thissnr=3.
        lowsnr=1.
        ovcheck=1
        while ovcheck>0:
            hemchoose=np.asarray([np.random.permutation(2)[0] for hemc in range(nseeds)])
            seedchoose=np.zeros(hemchoose.shape)
            for hemcnt, itrhem in enumerate(hemchoose):
                lhem=-1
                lhemperm=np.random.permutation(fwd_idxl)
                rhemperm=np.random.permutation(fwd_idxr)
                rhem=-1
                if itrhem==0:
                    lhem=lhem+1
                    seedchoose[hemcnt]=int(lhemperm[lhem])
                else:
                    rhem=rhem+1
                    seedchoose[hemcnt]=int(rhemperm[rhem])
            seedss=seedchoose.astype(int)
            print ovcheck
            seedov=np.zeros((len(parc_name),))
            labellist=mne.grow_labels(subject=subject_avg,seeds=seedss,extents=ext,hemis=hemchoose,subjects_dir=data_path)
            for which_parc in range(len(parc_name)):                
                labels_prec = mne.read_labels_from_annot('fsaverage', parc=parc_name[which_parc], subjects_dir=subjects_dir)
                if which_parc==0:
                    labels_prec.remove(labels_prec[-1])
                if which_parc==2:
                    labels_prec.remove(labels_prec[-1])
                    labels_prec.remove(labels_prec[-1])
                #labels_pre.remove(labels_pre[-1])
                #    labels_pre.remove(labels_pre[-1])
                
                path_Rmat=results_path+Rmat_names[which_parc]
                ltc=scipy.io.loadmat(path_Rmat,mat_dtype=True); ctf_to_ypos=np.squeeze(ltc['ctf_to_ypos'])            
                
                path_ypos=results_path+ypos_names[which_parc]
                ltc=scipy.io.loadmat(path_ypos,mat_dtype=True); ypos_sort=np.squeeze(ltc['ypos_sort'])
                labelsc=[labels_prec[jj] for jj in ypos_sort]
                for lbl in labelsc:
                    lbl.values.fill(1.0)
                labellist_sc=copy.deepcopy(labellist)
                labov_matc=np.zeros((len(labellist_sc),len(labelsc)))
                labels_mod=copy.deepcopy(labelsc)
                lpcnt=-1
                for lscnt,seedlab in enumerate(labellist_sc):
                    for lpcnt,parclab in enumerate(labelsc):
                        if seedlab.name[-2:]==parclab.name[-2:]:
                            labov_matc[lscnt,lpcnt]= len(np.intersect1d(seedlab.vertices, parclab.vertices))
    #                    else:
    #                        labov_mat[lscnt,lpcnt]=0                  
                iarprec=np.argmax(labov_matc,1)
                imc=np.max(labov_matc,1)
                iarc=iarprec.copy()
                iarposc=iarprec[imc>0]
                print len(np.unique(iarposc)),len(iarposc)
                print iarposc
                if len(np.unique(iarposc))<len(iarposc):
                   seedov[which_parc]=1 
                else:
                   seedov[which_parc]=0
            print seedov
            ovcheck=np.sum(seedov)
                
                
        parc0_connmat=np.zeros((parc_len[0],parc_len[0],10,len(subject_inds)))#nconxnconxnmeths,nsubs
        parc1_connmat=np.zeros((parc_len[1],parc_len[1],10,len(subject_inds)))
        parc2_connmat=np.zeros((parc_len[2],parc_len[2],10,len(subject_inds)))
        parc3_connmat=np.zeros((parc_len[3],parc_len[3],10,len(subject_inds)))
        parc4_connmat=np.zeros((parc_len[4],parc_len[4],10,len(subject_inds)))
        n_epochs=40
        #    LO_phase=np.pi*(2*np.random.random(n_epochs,)-1)
        #    RMF_phase=np.pi*(2*np.random.random(n_epochs,)-1)
        ROIs_phase=np.pi*(2*np.random.random((nseeds,n_epochs))-1)
        const_phase=np.pi*(2*np.random.random((nseeds,))-1)
        allfs=np.linspace(2,8,nseeds)
        mixfs=allfs[np.random.permutation(len(allfs))]
        fcoeff=mixfs[:nseeds]
        
        blen=125
        
        nconn=int(math.ceil(connperc*nseeds*(nseeds-1)/2))
        
        comseed=[1,2]
        comc=-1
        while len(comseed)>0 and comc<1000:
            comc=comc+1
            print comc        
            comseed=list()        
            Amat=np.eye(nseeds)#np.zeros((nseeds,nseeds))
            for ccnt in range(nconn):
                thisi=np.random.permutation(nseeds)[0]
                thisj=np.random.permutation(nseeds)[0]
                while Amat[thisi,thisj]==1:
                    thisi=np.random.permutation(nseeds)[0]
                    thisj=np.random.permutation(nseeds)[0]
                while thisi==thisj:
                    thisj=np.random.permutation(nseeds)[0]
                Amat[thisi,thisj]=1; Amat[thisj,thisi]=1;
            Amat=np.tril(Amat)
            for amatc1 in range(Amat.shape[0]): #rows
                for amatc2 in range(Amat.shape[1]):
                    if amatc1>amatc2 and Amat[amatc1,amatc2]==0:
                        rowovs=np.where(np.logical_and(Amat[amatc1,:]==Amat[amatc2,:],Amat[amatc2,:]==1))[0]
                        if len(rowovs)>0:
                            comseed.append(np.hstack((amatc1,amatc2,rowovs)))
                
            if comc==0:
                best_comseed=copy.deepcopy(comseed)
                best_Amat=Amat.copy()
            if len(comseed)<len(best_comseed):
                best_comseed=copy.deepcopy(comseed)
                best_Amat=Amat.copy()
        csl=[len(csl1) for csl1 in best_comseed]
        reqconn=nconn+nseeds
        extracons=np.sum(csl)-2*len(csl)
        actconn=reqconn+extracons
        comseed_mod=copy.deepcopy(best_comseed)
        comc=-1
        while len(comseed_mod)>0 and actconn>reqconn:# or np.min(csl)>3:
            comc=comc+1
            print comc
            this_targ=comseed_mod[np.argmin(csl)]
            thisizero=this_targ[np.random.permutation(2)[0]]
            thisjzero=this_targ[np.random.permutation(len(this_targ)-2)[0]+2]
            best_Amat[np.max(np.array([thisizero,thisjzero])),np.min(np.array([thisizero,thisjzero]))]=0.0
            comseed_mod=list()
            for amatc1 in range(best_Amat.shape[0]): #rows
                for amatc2 in range(best_Amat.shape[1]):
                    if amatc1>amatc2 and (best_Amat[amatc1,amatc2])==0:
                        rowovs=np.where(np.logical_and(best_Amat[amatc1,:]==best_Amat[amatc2,:],best_Amat[amatc2,:]==1))[0]
                        if len(rowovs)>0:
                            comseed_mod.append(np.hstack((amatc1,amatc2,rowovs)))
            csl=[len(csl1) for csl1 in comseed_mod]
            extracons=np.sum(csl)-2*len(csl)
            actconn=np.sum(best_Amat)+extracons
        Amat=best_Amat.copy()
        comseed=copy.deepcopy(comseed_mod)
        print Amat            
        t=np.linspace(-4*np.pi,4*np.pi,3020/5)
        datamat_grand=np.zeros((nseeds,blen+len(t),n_epochs))
        lowsnr_datamat_grand=np.zeros((nseeds,blen+len(t),n_epochs))
        datamat1=np.zeros((nseeds,blen+len(t)))
        lowsnr_datamat1=np.zeros((nseeds,blen+len(t)))
#        lowsnr_Amat=np.eye(Amat.shape[0],Amat.shape[1])
        for ecnt in range(n_epochs):
            print ecnt
            for rcnt in range(nseeds):
                seed_sig_forR=np.zeros((blen+len(t),nseeds))
                for scnt in range(nseeds):
                    seed_sig_forR[:,scnt]=np.concatenate((np.zeros(125,),1e-9*np.sin(fcoeff[scnt]*t+ROIs_phase[scnt,ecnt]+const_phase[rcnt])),axis=0)                  
                datamat1[rcnt,:]=np.matmul(Amat[rcnt,:],seed_sig_forR.transpose())#.transpose()  
                lowsnr_datamat1[rcnt,:]=datamat1[rcnt,:].copy()#np.matmul(nullAmat[rcnt,:],seed_sig_forR.transpose())
            datamat_grand[:,:,ecnt]=datamat1.copy()
            lowsnr_datamat_grand[:,:,ecnt]=lowsnr_datamat1.copy()
        ii=-1   
        actsrc_grandmat=list() 
        misscatch_grandmat=list() 
        actsrcfinal_grandmat=list()
#        idealconn_grandmat=list()
        missedcons_grandmat=list()
        actsrc_ovs=list()
        Amat_grandmat=list()
        list_allc=copy.deepcopy(list_all)
        list_incl=[list_allc[subcnt] for subcnt in subject_inds]
        for cnt, meg in enumerate(list_incl):
            ii=ii+1
            print cnt
            subjects_dir = data_path 
            fname_fwd = data_path + meg + 'ico5_forward_5-3L-EMEG-fwd.fif'
            fname_inv = inv_path + meg + 'InvOp_ico5newreg_fft_1_48_clean_ica_EMEG-inv.fif'#'InvOp_ico5diagreg_fft_1_48_clean_ica_EMEG-inv.fif'
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
                        
            src_sub=forward['src']
            leadfield=forward['sol']['data']
            
            for which_parc in range(len(parc_name)):#np.array([2,3]):#
                
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
                labels=[labels_pre[jj] for jj in ypos_sort]
                for lbl in labels:
                    lbl.values.fill(1.0)
                labellist_s=copy.deepcopy(labellist)
                labov_mat=np.zeros((len(labellist_s),len(labels)))
                labels_mod=copy.deepcopy(labels)
                lpcnt=-1
                for lscnt,seedlab in enumerate(labellist_s):
                    for lpcnt,parclab in enumerate(labels):
                        if seedlab.name[-2:]==parclab.name[-2:]:
                            labov_mat[lscnt,lpcnt]= len(np.intersect1d(seedlab.vertices, parclab.vertices))
    #                    else:
    #                        labov_mat[lscnt,lpcnt]=0                  
                iarpre=np.argmax(labov_mat,1)
                im=np.max(labov_mat,1)
                iar=iarpre.copy()
                iarpos=iarpre[im>0]
                iar[im==0]=1000
                misscatch=iarpre.copy()
                misscatch[im==0]=0
                misscatch[im>0]=1
                print iar
                missed_src=len(iarpre)-len(iar)
                actsrc_miss=iar.copy()
                
#                labellist_m=copy.deepcopy(labels)
#                morphed_labelsall_pre=[labellist_m[jj].morph(subject_from='fsaverage', smooth=5, subject_to=subjects[subject_no], subjects_dir=data_path) for jj in range(len(labellist_m))]
#                colors=nc(int(len(morphed_labelsall_pre)/2),bytes_=True, cmap='gist_rainbow')
#                colors=colors[np.random.permutation(int(len(morphed_labelsall_pre)/2))]
#                colors_split=colors[:int(len(morphed_labelsall_pre)/2),:]
                fname_split=parc_name[which_parc]+'_extsims'#+str(perm_inds[0])
#                make_parcellation_mine(src_sub,morphed_labelsall_pre,colors_split,subjects[subject_no],fname_split)
                ## write split labels
                
                morphed_labelsall = mne.read_labels_from_annot(subject=subjects[subject_no], parc=subjects[subject_no]+'_'+fname_split, subjects_dir=data_path)
#                fnamel=data_path + subjects[subject_no]+'/label/lh.'+subjects[subject_no]+'_'+fname_split+'.annot'    
#                fnamer=data_path + subjects[subject_no]+'/label/rh.'+subjects[subject_no]+'_'+fname_split+'.annot'
#                os.remove(fnamel)
#                os.remove(fnamer)
                print "hi"
                for lblcnt,labelcheck in enumerate(morphed_labelsall):
                    if labelcheck.hemi=='lh':
                        lblcs=np.intersect1d(src_sub[0]['vertno'],labelcheck.vertices)
                    else:
                        lblcs=np.intersect1d(src_sub[1]['vertno'],labelcheck.vertices)
                    if len(lblcs)<1:
                        print labelcheck
                        labellist_mc=copy.deepcopy(labels)
                        morphed_labelsall[lblcnt]=labellist_mc[lblcnt].morph(subject_from='fsaverage', smooth=5, subject_to=subjects[subject_no], subjects_dir=data_path)
    
                labelsp=copy.deepcopy(morphed_labelsall)
                pick_seeds=[labelsp[pscnt] for pscnt in iarpos]#fsaverage space
    #            morphed_labels=copy.deepcopy(pick_seeds)
                for thisps in pick_seeds:
                    thisps.values.fill(1.0)
                pick_seeds_cmas=np.asarray([llcm.center_of_mass(subject=subjects[subject_no],restrict_vertices=True, subjects_dir=data_path) for llcm in pick_seeds])
                pick_seeds_heml=list()
                for pscnt,llcm in enumerate(pick_seeds):
                    if llcm.hemi=='lh':
                        pick_seeds_heml.append(0)
                    else:
                        pick_seeds_heml.append(1)
                pick_seeds_hem=np.asarray(pick_seeds_heml)
                final_ext=5
                morphed_labels=mne.grow_labels(subject=subjects[subject_no],seeds=pick_seeds_cmas,extents=final_ext,hemis=pick_seeds_hem,subjects_dir=data_path,overlap=False)
    
    #            morphed_labels=[picksds[jj].morph(subject_from='fsaverage', smooth=5, subject_to=subjects[subject_no], subjects_dir=data_path) for jj in range(len(picksds))]
                labov_mat_final=np.zeros((len(morphed_labels),len(morphed_labelsall)))
                lpcnt1=-1
                for lscnt,seedlab in enumerate(morphed_labels):
                    for lpcnt1,parclab in enumerate(morphed_labelsall):
                        if seedlab.name[-2:]==parclab.name[-2:]:
                            labov_mat_final[lscnt,lpcnt1]= len(np.intersect1d(seedlab.vertices, parclab.vertices))
                labov_mat_final[labov_mat_final>0]=1
                ov_sum=np.sum(labov_mat_final,1)
                ov_sum[ov_sum<1]=1
                howmany_actsrc=np.sum(ov_sum)
                actsrc_final=np.where(labov_mat_final==1)[1]
                
                #labov_mat[labov_mat>0]=1
#                idealconn=np.eye(len(morphed_labelsall))
#                amatt=Amat.copy()
#                amatt=amatt[misscatch>0,:]
#                amatt=amatt[:,misscatch>0]
#                if len(iarpos)>1:
#                    for ccc1,iarc1 in enumerate(iarpos):
#                        for ccc2,iarc2 in enumerate(iarpos):
#                            if amatt[ccc1,ccc2]>0:
#                                idealconn[iarc1,iarc2]=1
#                for lovc in range(labov_mat_final.shape[0]):
#                    this_row=np.where(labov_mat_final[lovc,:]>0)[0]
#                    if len(this_row)>1:
#                        for ccc1,iarc1 in enumerate(this_row):
#                            for ccc2,iarc2 in enumerate(this_row):
#                                idealconn[iarc1,iarc2]=1
#                for lovc in range(labov_mat_final.shape[0]):
#                    this_row=np.where(labov_mat_final[lovc,:]>0)[0]
#                    for this_roi1 in this_row:
#                        for this_roi2 in this_row:
#                            for this_parcel in range(idealconn.shape[0]):
#                                if idealconn[this_roi1,this_parcel]==1:# or idealconn[this_parcel,this_roi1]>0:
#                                    idealconn[this_roi2,this_parcel]=1
#                                    idealconn[this_parcel,this_roi2]=1
                extramisses=0
                for seed_pairs in comseed:
#                    idealconn[seed_pairs[0],seed_pairs[1]]=1
#                    idealconn[seed_pairs[1],seed_pairs[0]]=1
                    if seed_pairs[0] in np.where(misscatch==0)[0] or seed_pairs[1] in np.where(misscatch==0)[0]:
                        extramisses=extramisses+1
                missedcons=np.tril(Amat)-np.eye(Amat.shape[0])
                #            missedcons=missedcons[misscatch==0,:]
                mct=missedcons[misscatch==0,:]
                mct=mct[:,misscatch==0]
                missedcons=np.sum(missedcons[misscatch==0,:])+np.sum(missedcons[:,misscatch==0])-np.sum(mct)+extramisses
                if cnt==0:
                    actsrc_grandmat.append(actsrc_miss)
                    actsrcfinal_grandmat.append(actsrc_final)
                    missedcons_grandmat.append(missedcons)
                    actsrc_ovs.append(ov_sum)
                    misscatch_grandmat.append(misscatch)
    #            labov_mat[labov_mat>0]=1
    #            ov_sum=np.sum(labov_mat,1)
    #            ov_sum[ov_sum<1]=1
    #            howmany_actsrc=np.sum(ov_sum)
    #            labellist_m=copy.deepcopy(pick_seeds)
    #[labellist_m[jj].morph(subject_from='fsaverage', smooth=5, subject_to=subjects[subject_no], subjects_dir=data_path) for jj in range(len(labellist_m))]
                parc_datamat_grand=datamat_grand[np.where(im>0)[0],:,:].copy()
                lowsnr_parc_datamat_grand=lowsnr_datamat_grand[np.where(im>0)[0],:,:].copy()
                label_colors = [label.color for label in labels]
                stc_sim_ideal=range(n_epochs)
                stc_sim_realideal=range(n_epochs)
                lowsnr_stc_sim=range(n_epochs)
                lowsnr_stc_sim_ideal=range(n_epochs)
                for ecnt in range(n_epochs):  
                    datamat=parc_datamat_grand[:,:,ecnt].copy()
                    lowsnr_datamat=lowsnr_parc_datamat_grand[:,:,ecnt].copy()
                    this_std=np.sqrt(np.mean(datamat[0,125:]**2)/(thissnr**2))#SNR=3 
                    lowsnr_this_std=np.sqrt(np.mean(lowsnr_datamat[0,125:]**2)/(lowsnr**2))
                    tmin=-0.124
                    tstep=0.001
                    lvert=src_sub[0]['vertno']
                    rvert=src_sub[1]['vertno']
                    vertices_sub = [lvert, rvert]
                    thislabel_data=np.zeros((len(lvert)+len(rvert),datamat.shape[1]))
#                    allnoise_thislabel_data=np.zeros((len(lvert)+len(rvert),datamat.shape[1]))
                    lowsnr_thislabel_data=np.zeros((len(lvert)+len(rvert),lowsnr_datamat.shape[1]))
               
                    for llcnt,ll1 in enumerate(morphed_labels):
                        if ll1.hemi == 'rh':
                            # for RH labels, add number of LH vertices
                            offset = forward['src'][0]['vertno'].shape[0]
                            # remember whether we are in the LH or RH
                            this_hemi = 1
                        elif ll1.hemi == 'lh':
                            offset = 0
                            this_hemi = 0
                
                        # get vertices on cortical surface inside label
                        idx = np.intersect1d(ll1.vertices, forward['src'][this_hemi]['vertno'])
                
                        # get vertices in source space inside label
                        fwd_idx = np.searchsorted(forward['src'][this_hemi]['vertno'], idx)
                        thislabel_data[fwd_idx+offset,:]=datamat[llcnt,:].copy()
                        lowsnr_thislabel_data[fwd_idx+offset,:]=lowsnr_datamat[llcnt,:].copy()
                    ## project from source space to sensor space
                    # multiply by leadfield
                    not_noise=np.unique(np.where(thislabel_data!=0)[0])
                    tobe_noise=np.asarray([thiszero for thiszero in range(thislabel_data.shape[0]) if thiszero not in not_noise])
#                    tobe_allnoise=np.asarray([thiszero for thiszero in range(thislabel_data.shape[0])])
    
                    noise_std=this_std.copy()#np.sqrt(np.mean(datamat[0,125:]**2))
                    lowsnr_noise_std=lowsnr_this_std.copy()#np.sqrt(np.mean(nulldatamat[0,125:]**2))
                    noise_mat=noise_std*np.random.randn(tobe_noise.shape[0],thislabel_data[:,125:].shape[1])+np.mean(datamat[0,125:])
                    lowsnr_noise_mat=lowsnr_noise_std*np.random.randn(tobe_noise.shape[0],lowsnr_thislabel_data[:,125:].shape[1])+np.mean(lowsnr_datamat[0,125:])
#                    allnoise_mat=noise_std*np.random.randn(tobe_allnoise.shape[0],thislabel_data[:,125:].shape[1])+np.mean(datamat[0,125:])
    
                    thislabel_data_ideal=thislabel_data.copy()
                    thislabel_data_ideal[tobe_noise,125:]=0.001*noise_mat
                    thislabel_data[tobe_noise,125:]=noise_mat
                    lowsnr_thislabel_data[tobe_noise,125:]=lowsnr_noise_mat
#                    allnoise_thislabel_data[:,125:]=allnoise_mat
                    print ecnt
                    rmat=this_std*np.random.randn(thislabel_data.shape[0],thislabel_data.shape[1])
                    lowsnr_rmat=lowsnr_this_std*np.random.randn(lowsnr_thislabel_data.shape[0],lowsnr_thislabel_data.shape[1])
    
                    thislabel_datan=thislabel_data+rmat
#                    allnoise_datan=allnoise_thislabel_data+rmat
                    lowsnr_thislabel_datan=lowsnr_thislabel_data+lowsnr_rmat
                    stc_sim_ideal[ecnt]=mne.SourceEstimate(thislabel_datan, vertices=vertices_sub, tmin=tmin, tstep=tstep, subject=subjects[subject_no], verbose=None)#thislabel_datan
                    stc_sim_realideal[ecnt]=mne.SourceEstimate(thislabel_data_ideal, vertices=vertices_sub, tmin=tmin, tstep=tstep, subject=subjects[subject_no], verbose=None)#thislabel_datan
                    lowsnr_stc_sim_ideal[ecnt]=mne.SourceEstimate(lowsnr_thislabel_datan, vertices=vertices_sub, tmin=tmin, tstep=tstep, subject=subjects[subject_no], verbose=None)#thislabel_datan
    
                    evoked_data=np.matmul(leadfield,thislabel_datan)[np.newaxis,:,:]
                    lowsnr_evoked_data=np.matmul(leadfield,lowsnr_thislabel_datan)[np.newaxis,:,:]
#                    allnoise_evoked_data=np.matmul(leadfield,allnoise_datan)[np.newaxis,:,:]
                    if ecnt==0:
                        epochs_data=evoked_data.copy()
                        lowsnr_epochs_data=lowsnr_evoked_data.copy()
#                        allnoise_epochs_data=allnoise_evoked_data.copy()
                        events=np.array([0,1,1])[np.newaxis,:]
            
                    else:
                        epochs_data=np.append(epochs_data,evoked_data,axis=0)
#                        allnoise_epochs_data=np.append(allnoise_epochs_data,allnoise_evoked_data,axis=0)
                        lowsnr_epochs_data=np.append(lowsnr_epochs_data,lowsnr_evoked_data,axis=0)
                        this_event=np.array([ecnt+1,1,1])[np.newaxis,:]
                        events=np.append(events,this_event,axis=0)
                # create evoked array
            #    evoked_fwd = mne.EvokedArray(evoked_data, info=info, tmin=-0.769)
                epochs_fwd = mne.EpochsArray(epochs_data, info=epochs1.info, events=events, event_id=1, tmin=tmin)
#                allnoise_epochs_fwd = mne.EpochsArray(allnoise_epochs_data, info=epochs1.info, events=events, event_id=1, tmin=tmin)
                lowsnr_epochs_fwd = mne.EpochsArray(lowsnr_epochs_data, info=epochs1.info, events=events, event_id=1, tmin=tmin)
    
                snr = thissnr#3.0#ediv.mean()
                lambda2 = 1.0 / snr ** 2
                lambda2lowsnr = 1.0 / lowsnr ** 2
                method = 'MNE'
                ## apply inverse to the evoked array
                #stc_sim = apply_inverse(evoked_fwd, inverse_operator_eegmeg, lambda2,method=method, pick_ori="normal")
                stc_sim = apply_inverse_epochs(epochs_fwd, inverse_operator_eegmeg, lambda2,method=method, pick_ori="normal")
#                allnoise_stc_sim = apply_inverse_epochs(allnoise_epochs_fwd, inverse_operator_eegmeg, lambda2,method=method, pick_ori="normal")
                lowsnr_stc_sim = apply_inverse_epochs(lowsnr_epochs_fwd, inverse_operator_eegmeg, lambda2lowsnr,method=method, pick_ori="normal")
                print "extracting label time course"
                label_ts = mne.extract_label_time_course(stc_sim, morphed_labelsall, src_sub, mode='mean_flip',return_generator=False)
#                allnoise_label_ts = mne.extract_label_time_course(allnoise_stc_sim, morphed_labelsall, src_sub, mode='mean_flip',return_generator=False)
                lowsnr_label_ts = mne.extract_label_time_course(lowsnr_stc_sim, morphed_labelsall, src_sub, mode='mean_flip',return_generator=False)
                lowsnr_label_ts_ideal = mne.extract_label_time_course(lowsnr_stc_sim_ideal, morphed_labelsall, src_sub, mode='mean_flip',return_generator=False)
                label_ts_ideal = mne.extract_label_time_course(stc_sim_ideal, morphed_labelsall, src_sub, mode='mean_flip',return_generator=False)
                label_ts_realideal = mne.extract_label_time_course(stc_sim_realideal, morphed_labelsall, src_sub, mode='mean_flip',return_generator=False)
        
                    ##################################
                    
                #    cl=range(len(labels))
                #    clw=range(len(labels))
                #    cl_ideal=range(len(labels))
                #    clw_ideal=range(len(labels))
                #    cl_realideal=range(len(labels))
                #    clw_realideal=range(len(labels))
                #    for cntlts in range(n_epochs):
                #        for cntcl in range(len(labels)):        
                #            cl[cntcl]=simstc_morphed[cntlts].in_label(labels[cntcl]).data
                #            clw[cntcl]=np.argmax(np.mean(np.abs(cl[cntcl][:,125:]),axis=1))
                #            thislabelsign=mne.label.label_sign_flip(labels[cntcl],src_avg)[clw[cntcl]]
                #            label_ts[cntlts][cntcl,:]=thislabelsign*np.expand_dims(cl[cntcl][clw[cntcl],:].transpose(),0)
                #            ##ideal
                #            cl_ideal[cntcl]=simstc_morphed_ideal[cntlts].in_label(labels[cntcl]).data
                #            clw_ideal[cntcl]=np.argmax(np.mean(np.abs(cl_ideal[cntcl][:,125:]),axis=1))
                #            thislabelsign2=mne.label.label_sign_flip(labels[cntcl],src_avg)[clw_ideal[cntcl]]
                #            label_ts_ideal[cntlts][cntcl,:]=thislabelsign2*np.expand_dims(cl_ideal[cntcl][clw_ideal[cntcl],:].transpose(),0)
                #            ##realideal
                #            cl_realideal[cntcl]=simstc_morphed_realideal[cntlts].in_label(labels[cntcl]).data
                #            clw_realideal[cntcl]=np.argmax(np.mean(np.abs(cl_realideal[cntcl][:,125:]),axis=1))
                #            thislabelsign2=mne.label.label_sign_flip(labels[cntcl],src_avg)[clw_realideal[cntcl]]
                #            label_ts_realideal[cntlts][cntcl,:]=thislabelsign2*np.expand_dims(cl_realideal[cntcl][clw_realideal[cntcl],:].transpose(),0)
        
                print "extracting label time course DONE"
            
                ##coherence
                
                fmin = 8.
                fmax = 55.
                sfreq = raw.info['sfreq']  # the sampling frequency
                con_methods = ['coh']#['pli2_unbiased']
                print "computing first coh"
                con, freqs, times, n_epochs, n_tapers = spectral_connectivity(
                label_ts, method=con_methods, mode='multitaper', sfreq=sfreq, fmin=fmin,fmax=fmax, 
                faverage=True, mt_adaptive=True, n_jobs=1)                
                con_res = np.squeeze(con.copy())#np.mean(con,axis=2)#dict()
                print "computing second coh"
                con_ideal, freqs, times, n_epochs, n_tapers = spectral_connectivity(
                    label_ts_ideal, method=con_methods, mode='multitaper', sfreq=sfreq, fmin=fmin,
                    fmax=fmax, faverage=True, mt_adaptive=True, n_jobs=1)
                con_res_ideal = np.squeeze(con_ideal.copy())#np.mean(con_ideal,axis=2)#con_ideal[:,:,0]#dict()
                print "computing 3rd coh"
                con_realideal, freqs, times, n_epochs, n_tapers = spectral_connectivity(
                    label_ts_realideal, method=con_methods, mode='multitaper', sfreq=sfreq, fmin=fmin,
                    fmax=fmax, faverage=True, mt_adaptive=True, n_jobs=1)
                con_res_realideal = np.squeeze(con_realideal.copy())
                print "computing 4th coh"
                lowsnr_con, freqs, times, n_epochs, n_tapers = spectral_connectivity(
                lowsnr_label_ts, method=con_methods, mode='multitaper', sfreq=sfreq, fmin=fmin,fmax=fmax, 
                faverage=True, mt_adaptive=True, n_jobs=1)                
                lowsnr_con_res = np.squeeze(lowsnr_con.copy())
                print "computing 5th coh"
                lowsnr_con_ideal, freqs, times, n_epochs, n_tapers = spectral_connectivity(
                lowsnr_label_ts_ideal, method=con_methods, mode='multitaper', sfreq=sfreq, fmin=fmin,fmax=fmax, 
                faverage=True, mt_adaptive=True, n_jobs=1)                
                lowsnr_con_res_ideal = np.squeeze(lowsnr_con_ideal.copy())#np.mean(con,axis=2)#dict()
                
                print "computing 1st imcoh"
                con_methods2 = ['imcoh']#['pli2_unbiased']
                con2,  freqs, times, n_epochs, n_tapers = spectral_connectivity(
                    label_ts, method=con_methods2, mode='multitaper', sfreq=sfreq, fmin=fmin,
                    fmax=fmax, faverage=True, mt_adaptive=True, n_jobs=1)
            
                con_res2 = np.squeeze(con2.copy())#np.mean(con2,axis=2)#con2[:,:,0]#dict()
                print "computing second imcoh"
                con_ideal2, freqs, times, n_epochs, n_tapers = spectral_connectivity(
                    label_ts_ideal, method=con_methods2, mode='multitaper', sfreq=sfreq, fmin=fmin,
                    fmax=fmax, faverage=True, mt_adaptive=True, n_jobs=1)
                #     con is a 3D array, get the connectivity for the first (and only) freq. band
                #     for each method
                con_res_ideal2 = np.squeeze(con_ideal2.copy())#np.mean(con_ideal2,axis=2)#con_ideal2[:,:,0]#dict()
                print "computing 3rd imcoh"
                con_realideal2, freqs, times, n_epochs, n_tapers = spectral_connectivity(
                    label_ts_realideal, method=con_methods2, mode='multitaper', sfreq=sfreq, fmin=fmin,
                    fmax=fmax, faverage=True, mt_adaptive=True, n_jobs=1)
                    
                print "computing 4th coh"
                lowsnr_con2, freqs, times, n_epochs, n_tapers = spectral_connectivity(
                lowsnr_label_ts, method=con_methods2, mode='multitaper', sfreq=sfreq, fmin=fmin,fmax=fmax, 
                faverage=True, mt_adaptive=True, n_jobs=1)                
                lowsnr_con_res2 = np.squeeze(lowsnr_con2.copy())
                print "computing 5th coh"
                lowsnr_con_ideal2, freqs, times, n_epochs, n_tapers = spectral_connectivity(
                lowsnr_label_ts_ideal, method=con_methods2, mode='multitaper', sfreq=sfreq, fmin=fmin,fmax=fmax, 
                faverage=True, mt_adaptive=True, n_jobs=1)                
                lowsnr_con_res_ideal2 = np.squeeze(lowsnr_con_ideal2.copy())
                #     con is a 3D array, get the connectivity for the first (and only) freq. band
                #     for each method
                con_res_realideal2 = np.squeeze(con_realideal2.copy())
                if which_parc==0:
                    parc0_connmat[:,:,0,cnt]=con_res.copy()#con_res_final=con_res[:,:,np.newaxis].copy()
                    parc0_connmat[:,:,1,cnt]=con_res_ideal.copy()
                    parc0_connmat[:,:,2,cnt]=con_res_realideal.copy()
                    parc0_connmat[:,:,3,cnt]=lowsnr_con_res.copy()
                    parc0_connmat[:,:,4,cnt]=lowsnr_con_res_ideal.copy()
                    parc0_connmat[:,:,5,cnt]=con_res2.copy()
                    parc0_connmat[:,:,6,cnt]=con_res_ideal2.copy()
                    parc0_connmat[:,:,7,cnt]=con_res_realideal2.copy()
                    parc0_connmat[:,:,8,cnt]=lowsnr_con_res2.copy()
                    parc0_connmat[:,:,9,cnt]=lowsnr_con_res_ideal2.copy()
                if which_parc==1:
                    parc1_connmat[:,:,0,cnt]=con_res.copy()#con_res_final=con_res[:,:,np.newaxis].copy()
                    parc1_connmat[:,:,1,cnt]=con_res_ideal.copy()
                    parc1_connmat[:,:,2,cnt]=con_res_realideal.copy()
                    parc1_connmat[:,:,3,cnt]=lowsnr_con_res.copy()
                    parc1_connmat[:,:,4,cnt]=lowsnr_con_res_ideal.copy()
                    parc1_connmat[:,:,5,cnt]=con_res2.copy()
                    parc1_connmat[:,:,6,cnt]=con_res_ideal2.copy()
                    parc1_connmat[:,:,7,cnt]=con_res_realideal2.copy()
                    parc1_connmat[:,:,8,cnt]=lowsnr_con_res2.copy()
                    parc1_connmat[:,:,9,cnt]=lowsnr_con_res_ideal2.copy()
                if which_parc==2:
                    parc2_connmat[:,:,0,cnt]=con_res.copy()#con_res_final=con_res[:,:,np.newaxis].copy()
                    parc2_connmat[:,:,1,cnt]=con_res_ideal.copy()
                    parc2_connmat[:,:,2,cnt]=con_res_realideal.copy()
                    parc2_connmat[:,:,3,cnt]=lowsnr_con_res.copy()
                    parc2_connmat[:,:,4,cnt]=lowsnr_con_res_ideal.copy()
                    parc2_connmat[:,:,5,cnt]=con_res2.copy()
                    parc2_connmat[:,:,6,cnt]=con_res_ideal2.copy()
                    parc2_connmat[:,:,7,cnt]=con_res_realideal2.copy()
                    parc2_connmat[:,:,8,cnt]=lowsnr_con_res2.copy()
                    parc2_connmat[:,:,9,cnt]=lowsnr_con_res_ideal2.copy()
                if which_parc==3:
                    parc3_connmat[:,:,0,cnt]=con_res.copy()#con_res_final=con_res[:,:,np.newaxis].copy()
                    parc3_connmat[:,:,1,cnt]=con_res_ideal.copy()
                    parc3_connmat[:,:,2,cnt]=con_res_realideal.copy()
                    parc3_connmat[:,:,3,cnt]=lowsnr_con_res.copy()
                    parc3_connmat[:,:,4,cnt]=lowsnr_con_res_ideal.copy()
                    parc3_connmat[:,:,5,cnt]=con_res2.copy()
                    parc3_connmat[:,:,6,cnt]=con_res_ideal2.copy()
                    parc3_connmat[:,:,7,cnt]=con_res_realideal2.copy()
                    parc3_connmat[:,:,8,cnt]=lowsnr_con_res2.copy()
                    parc3_connmat[:,:,9,cnt]=lowsnr_con_res_ideal2.copy()
                if which_parc==4:
                    parc4_connmat[:,:,0,cnt]=con_res.copy()#con_res_final=con_res[:,:,np.newaxis].copy()
                    parc4_connmat[:,:,1,cnt]=con_res_ideal.copy()
                    parc4_connmat[:,:,2,cnt]=con_res_realideal.copy()
                    parc4_connmat[:,:,3,cnt]=lowsnr_con_res.copy()
                    parc4_connmat[:,:,4,cnt]=lowsnr_con_res_ideal.copy()
                    parc4_connmat[:,:,5,cnt]=con_res2.copy()
                    parc4_connmat[:,:,6,cnt]=con_res_ideal2.copy()
                    parc4_connmat[:,:,7,cnt]=con_res_realideal2.copy()
                    parc4_connmat[:,:,8,cnt]=lowsnr_con_res2.copy()
                    parc4_connmat[:,:,9,cnt]=lowsnr_con_res_ideal2.copy()
        
        this_pathto=out_path+outname[0]+'parc_connmat_sim'+str(perm_inds[0]+sim_offset)+'.mat'#'DKA68_orig_modified_coh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot            
        scipy.io.savemat(this_pathto,{'parc0_connmat':parc0_connmat})
        this_pathto=out_path+outname[1]+'parc_connmat_sim'+str(perm_inds[0]+sim_offset)+'.mat'#'DKA68_orig_modified_coh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot            
        scipy.io.savemat(this_pathto,{'parc1_connmat':parc1_connmat})
        this_pathto=out_path+outname[2]+'parc_connmat_sim'+str(perm_inds[0]+sim_offset)+'.mat'#'DKA68_orig_modified_coh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot            
        scipy.io.savemat(this_pathto,{'parc2_connmat':parc2_connmat})
        this_pathto=out_path+outname[3]+'parc_connmat_sim'+str(perm_inds[0]+sim_offset)+'.mat'#'DKA68_orig_modified_coh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot            
        scipy.io.savemat(this_pathto,{'parc3_connmat':parc3_connmat})
        this_pathto=out_path+outname[4]+'parc_connmat_sim'+str(perm_inds[0]+sim_offset)+'.mat'#'DKA68_orig_modified_coh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot            
        scipy.io.savemat(this_pathto,{'parc4_connmat':parc4_connmat})
        import pickle
        this_pathto=out_path+'/activesources_grandmat_sim'+str(perm_inds[0]+sim_offset)+'.py'#'DKA68_orig_modified_coh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot            
        with open(this_pathto,'wb') as fp:
            pickle.dump(actsrc_grandmat,fp)
        this_pathto=out_path+'/activesources_howmanyovs_sim'+str(perm_inds[0]+sim_offset)+'.py'#'DKA68_orig_modified_coh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot            
        with open(this_pathto,'wb') as fp:
            pickle.dump(actsrc_ovs,fp)        
        this_pathto=out_path+'/activesourcesfinal_grandmat_sim'+str(perm_inds[0]+sim_offset)+'.py'#'DKA68_orig_modified_coh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot            
        with open(this_pathto,'wb') as fp:
            pickle.dump(actsrcfinal_grandmat,fp)
        this_pathto=out_path+'/missedcons_grandmat_sim'+str(perm_inds[0]+sim_offset)+'.py'#'DKA68_orig_modified_coh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot            
        with open(this_pathto,'wb') as fp:
            pickle.dump(missedcons_grandmat,fp)
        this_pathto=out_path+'/misscatch_grandmat_sim'+str(perm_inds[0]+sim_offset)+'.py'#'DKA68_orig_modified_coh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot            
        with open(this_pathto,'wb') as fp:
            pickle.dump(misscatch_grandmat,fp)
    
#with open ('outfile', 'rb') as fp:
#    itemlist = pickle.load(fp)
# First, we reorder the labels based on their location in the left hemi
#tfin=time.time()-t0
#label_names = [label.name for label in labels]
#
#lh_labels = [name for name in label_names if name.endswith('lh')]
#
## Get the y-location of the label
#label_ypos = list()
#for name in lh_labels:
#    idx = label_names.index(name)
#    ypos = np.mean(labels[idx].pos[:, 1])
#    label_ypos.append(ypos)
#
## Reorder the labels based on their location
#lh_labels = [label for (yp, label) in sorted(zip(label_ypos, lh_labels))]
#
## For the right hemi
#rh_labels = [label[:-2] + 'rh' for label in lh_labels]
#
## Save the plot order and create a circular layout
#node_order = list()
#node_order.extend(lh_labels[::-1])  # reverse the order
#node_order.extend(rh_labels)
#
#node_angles = circular_layout(label_names, node_order, start_pos=90,
#                              group_boundaries=[0, len(label_names) / 2])
#
## Plot the graph using node colors from the FreeSurfer parcellation. We only
## show the 300 strongest connections.
#plt.close("all")
#con_res_mean=np.squeeze(np.mean(con_res_final,axis=2))
#plot_connectivity_circle(np.abs(con_res_mean), range(len(label_names)),n_lines=None,
#                         node_angles=node_angles, node_colors=label_colors,#vmin=0,vmax=1,
#                         title='All-to-All Connectivity '
#                               'Coh')
##con_mod_mean=np.mean(con_mod_final,axis=2)                             
##plot_connectivity_circle(np.abs(con_mod_mean), range(len(label_names)),n_lines=300,
##                         node_angles=node_angles, node_colors=label_colors,#vmin=0,vmax=1e-19,
##                         title='All-to-All Connectivity '
##                               'Coh')
#con_res_mean_ideal=np.squeeze(np.mean(con_res_final_ideal,axis=2))
#plot_connectivity_circle(np.abs(con_res_mean_ideal), label_names, n_lines=None,#vmin=0,#vmax=1,
#                         node_angles=node_angles, node_colors=label_colors,
#                         title='All-to-All Connectivity no leakage Coh')
####
#con_res_mean_realideal=np.squeeze(np.mean(con_res_final_realideal,axis=2))
#plot_connectivity_circle(np.abs(con_res_mean_realideal), label_names, n_lines=None,#vmin=0,#vmax=1,
#                         node_angles=node_angles, node_colors=label_colors,
#                         title='All-to-All Connectivity '
#                               'Ideal Coh')                               
####
#con_res_mean2=np.squeeze(np.mean(con_res_final2,axis=2))
#plot_connectivity_circle(np.abs(con_res_mean2), label_names, n_lines=None,
#                         node_angles=node_angles, node_colors=label_colors,#vmin=0,vmax=1,
#                         title='All-to-All Connectivity '
#                               'ImCoh')
##con2_mod_mean=np.mean(con2_mod_final,axis=2)                             
##plot_connectivity_circle(np.abs(con2_mod_mean), range(len(label_names)),n_lines=300,
##                         node_angles=node_angles, node_colors=label_colors,#vmin=0,vmax=1,
##                         title='All-to-All Connectivity '
##                               'Coh')
#con_res_mean_ideal2=np.squeeze(np.mean(con_res_final_ideal2,axis=2))
#plot_connectivity_circle(np.abs(con_res_mean_ideal2), label_names, n_lines=None,#vmin=0,vmax=1,
#                         node_angles=node_angles, node_colors=label_colors,
#                         title='All-to-All Connectivity '
#                               'no leakage ImCoh') 
#                               
####
#con_res_mean_realideal2=np.squeeze(np.mean(con_res_final_realideal2,axis=2))
#plot_connectivity_circle(np.abs(con_res_mean_realideal2), label_names, n_lines=None,#vmin=0,vmax=1,
#                         node_angles=node_angles, node_colors=label_colors,
#                         title='All-to-All Connectivity '
#                               'Ideal ImCoh')                                
##plt.savefig('circle.png', facecolor='black')
#outname=['aparc','aparc_mod','aparc2009','aparc2009_mod','RG']
##out_path='/imaging/rf02/parcellation/extended_sim/'
##path_to=label_path+outname[which_parc]+'_coh_test.mat'#'DKA68_orig_modified_coh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot
##scipy.io.savemat(path_to,{'con_res_final':con_res_final})
###path_to=label_path+outname[which_parc]+'_coh_mod_test.mat'#'DKA68_orig_modified_coh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot
###scipy.io.savemat(path_to,{'con_mod_final':con_mod_final})
##
##path_to=label_path+outname[which_parc]+'_ideal_coh_test.mat'#'DKA68_ideal_modified_coh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot
##scipy.io.savemat(path_to,{'con_res_final_ideal':con_res_final_ideal})
##
##path_to=label_path+outname[which_parc]+'_imcoh_test.mat'#'DKA68_orig_modified_imcoh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot
##scipy.io.savemat(path_to,{'con_res_final2':con_res_final2})
###path_to=label_path+outname[which_parc]+'_imcoh_mod_test.mat'#'DKA68_orig_modified_imcoh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot
###scipy.io.savemat(path_to,{'con2_mod_final':con2_mod_final})
##
##path_to=label_path+outname[which_parc]+'_ideal_imcoh_test.mat'#'DKA68_ideal_modified_imcoh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot
##scipy.io.savemat(path_to,{'con_res_final_ideal2':con_res_final_ideal2})
#plt.show()
# 
#def mean2(x):
#    y=np.sum(x)/np.size(x)
#    return y
#def corr2(a,b):
#    a=a-mean2(a)
#    b=b-mean2(b)
#    r=(a*b).sum()/np.sqrt((a*a).sum()*(b*b).sum())
#    return r       
#cc2=(con_res_mean_realideal+con_res_mean_realideal.T)/2
##cc2=cc2+np.eye(cc2.shape[0])
#cc1=(con_res_mean+con_res_mean.T)/2 
##cc1=cc1+np.eye(cc1.shape[0])  
#R=corr2(cc1,cc2)