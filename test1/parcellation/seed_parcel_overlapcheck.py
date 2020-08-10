# -*- coding: utf-8 -*-
"""
Created on Sat May 20 17:23:55 2017

@author: rf02
"""

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



ll = range(10)#[0]
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
ext=5
connperc=0.25


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
    

labels_glasser = mne.read_labels_from_annot('fsaverage', parc='HCP-MMP1', subjects_dir=subjects_dir)
labels_glasser.remove(labels_glasser[0])
labels_glasser.remove(labels_glasser[0])            
            
            
            
#labels_pre.remove(labels_pre[-1])
#    labels_pre.remove(labels_pre[-1])
actsrc_grandmat=np.zeros((5,10,len(ll)))
for simcnt in ll:
    for nseeds in np.array([10]): #,10,15,20
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
#            hemchoose=np.asarray([np.random.permutation(2)[0] for hemc in range(nseeds)])
#            seedchoose=np.zeros(hemchoose.shape)
#            for hemcnt, itrhem in enumerate(hemchoose):
#                lhem=-1
#                lhemperm=np.random.permutation(fwd_idxl)
#                rhemperm=np.random.permutation(fwd_idxr)
#                rhem=-1
#                if itrhem==0:
#                    lhem=lhem+1
#                    seedchoose[hemcnt]=int(lhemperm[lhem])
#                else:
#                    rhem=rhem+1
#                    seedchoose[hemcnt]=int(rhemperm[rhem])
#            seedss=seedchoose.astype(int)
            print ovcheck
            seedov=np.zeros((len(parc_name),))
#            labellist=mne.grow_labels(subject=subject_avg,seeds=seedss,extents=ext,hemis=hemchoose,subjects_dir=data_path)
            glasserperm=np.random.permutation(len(labels_glasser))[:nseeds]
            lgc=copy.deepcopy(labels_glasser)
            labellist=[lgc[lgcn] for lgcn in range(len(lgc)) if lgcn in glasserperm]
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
#        actsrc_grandmat=list() 
#        misscatch_grandmat=list() 
#        actsrcfinal_grandmat=list()
##        idealconn_grandmat=list()
#        missedcons_grandmat=list()
        actsrc_ovs=list()
        Amat_grandmat=list()
        list_allc=copy.deepcopy(list_all)
        list_incl=[list_allc[subcnt] for subcnt in subject_inds]
        ii=-1
        
            ## simulate activity in labels
                        
       
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
            actsrc_grandmat[which_parc,:,simcnt]=actsrc_miss.copy()
agpre=actsrc_grandmat.copy()
actsrc_grandmat[actsrc_grandmat<1000]=0
actsrc_grandmat[actsrc_grandmat>0]=1
mm=np.mean(actsrc_grandmat,axis=1)
mmm=np.mean(mm,1)
#            actsrcfinal_grandmat.append(actsrc_final)
#            missedcons_grandmat.append(missedcons)
#            actsrc_ovs.append(ov_sum)
#            misscatch_grandmat.append(misscatch)               

#                fname_split=parc_name[which_parc]+'_extsims'#+str(perm_inds[0])
#                
#                morphed_labelsall = mne.read_labels_from_annot(subject=subjects[subject_no], parc=subjects[subject_no]+'_'+fname_split, subjects_dir=data_path)
#
#                print "hi"
#                for lblcnt,labelcheck in enumerate(morphed_labelsall):
#                    if labelcheck.hemi=='lh':
#                        lblcs=np.intersect1d(src_sub[0]['vertno'],labelcheck.vertices)
#                    else:
#                        lblcs=np.intersect1d(src_sub[1]['vertno'],labelcheck.vertices)
#                    if len(lblcs)<1:
#                        print labelcheck
#                        labellist_mc=copy.deepcopy(labels)
#                        morphed_labelsall[lblcnt]=labellist_mc[lblcnt].morph(subject_from='fsaverage', smooth=5, subject_to=subjects[subject_no], subjects_dir=data_path)
    
#                labelsp=copy.deepcopy(morphed_labelsall)
#                pick_seeds=[labelsp[pscnt] for pscnt in iarpos]#fsaverage space
#    #            morphed_labels=copy.deepcopy(pick_seeds)
#                for thisps in pick_seeds:
#                    thisps.values.fill(1.0)
#                pick_seeds_cmas=np.asarray([llcm.center_of_mass(subject=subjects[subject_no],restrict_vertices=True, subjects_dir=data_path) for llcm in pick_seeds])
#                pick_seeds_heml=list()
#                for pscnt,llcm in enumerate(pick_seeds):
#                    if llcm.hemi=='lh':
#                        pick_seeds_heml.append(0)
#                    else:
#                        pick_seeds_heml.append(1)
#                pick_seeds_hem=np.asarray(pick_seeds_heml)
#                final_ext=5
#                morphed_labels=mne.grow_labels(subject=subjects[subject_no],seeds=pick_seeds_cmas,extents=final_ext,hemis=pick_seeds_hem,subjects_dir=data_path,overlap=False)
    
    #            morphed_labels=[picksds[jj].morph(subject_from='fsaverage', smooth=5, subject_to=subjects[subject_no], subjects_dir=data_path) for jj in range(len(picksds))]
#                morphed_labels=copy.deepcopy(labellist)
#                morphed_labelsall=copy.deepcopy(labels)
#                labov_mat_final=np.zeros((len(morphed_labels),len(morphed_labelsall)))
#                lpcnt1=-1
#                for lscnt,seedlab in enumerate(morphed_labels):
#                    for lpcnt1,parclab in enumerate(morphed_labelsall):
#                        if seedlab.name[-2:]==parclab.name[-2:]:
#                            labov_mat_final[lscnt,lpcnt1]= len(np.intersect1d(seedlab.vertices, parclab.vertices))
#                labov_mat_final[labov_mat_final>0]=1
#                ov_sum=np.sum(labov_mat_final,1)
#                ov_sum[ov_sum<1]=1
#                howmany_actsrc=np.sum(ov_sum)
#                actsrc_final=np.where(labov_mat_final==1)[1]
#
#                extramisses=0
#                for seed_pairs in comseed:
##                    idealconn[seed_pairs[0],seed_pairs[1]]=1
##                    idealconn[seed_pairs[1],seed_pairs[0]]=1
#                    if seed_pairs[0] in np.where(misscatch==0)[0] or seed_pairs[1] in np.where(misscatch==0)[0]:
#                        extramisses=extramisses+1
#                missedcons=np.tril(Amat)-np.eye(Amat.shape[0])
#                #            missedcons=missedcons[misscatch==0,:]
#                mct=missedcons[misscatch==0,:]
#                mct=mct[:,misscatch==0]
#                missedcons=np.sum(missedcons[misscatch==0,:])+np.sum(missedcons[:,misscatch==0])-np.sum(mct)+extramisses
                
