# -*- coding: utf-8 -*-
"""
Created on Sun Oct  1 11:28:20 2017

@author: rf02
"""

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


sys.path.insert(1,'/imaging/rf02/scikit-learn-0.15.0')
sys.path.insert(1,'/home/rf02/.local/lib/python2.7/site-packages')
import sklearn
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
sys.path.insert(1,'/imaging/local/software/python_packages/scikit-learn/v0.13.1/lib.linux-x86_64-2.7/sklearn/externals')
import joblib
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
import scipy
from scipy.spatial.distance import pdist, squareform
from scipy import sparse
from scipy import stats as stats

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
import os
from mne.viz import circular_layout, plot_connectivity_circle


def get_components(A, no_depend=False):
    '''
    Returns the components of an undirected graph specified by the binary and
    undirected adjacency matrix adj. Components and their constitutent nodes
    are assigned the same index and stored in the vector, comps. The vector,
    comp_sizes, contains the number of nodes beloning to each component.
    Parameters
    ----------
    A : NxN np.ndarray
    binary undirected adjacency matrix
    no_depend : Any
    Does nothing, included for backwards compatibility
    Returns
    -------
    comps : Nx1 np.ndarray
    vector of component assignments for each node
    comp_sizes : Mx1 np.ndarray
    vector of component sizes
    Notes
    -----
    Note: disconnected nodes will appear as components with a component
    size of 1
    Note: The identity of each component (i.e. its numerical value in the
    result) is not guaranteed to be identical the value returned in BCT,
    matlab code, although the component topology is.
    Many thanks to Nick Cullen for providing this implementation
    '''
#    if not np.all(A == A.T): # ensure matrix is undirected
#        raise BCTParamError('get_components can only be computed for undirected'
#        ' matrices. If your matrix is noisy, correct it with np.around')
#    A = binarize(A, copy=True)
    n = len(A)
    np.fill_diagonal(A, 1)
    edge_map = [{u,v} for u in range(n) for v in range(n) if A[u,v] == 1]
    union_sets = []
    for item in edge_map:
        temp = []
        for s in union_sets:
            if not s.isdisjoint(item):
                item = s.union(item)
            else:
                temp.append(s)
        temp.append(item)
        union_sets = temp
    comps = np.array([i+1 for v in range(n) for i in range(len(union_sets)) if v in union_sets[i]])
    comp_sizes = np.array([len(s) for s in union_sets])
    return comps, comp_sizes


data_path = '/imaging/rf02/Semnet/' # root directory for your MEG data
os.chdir(data_path)
event_path = '/imaging/rf02/Semnet/'#'/group/erp/data/olaf.hauk/MEG/TypLexM/data/MF22_new/'
#label_path = '/group/erp/data/olaf.hauk/MEG/TypLexM/labels/createdlabels_SL/'
inv_path = '/imaging/rf02/Semnet/'
out_path = '/imaging/rf02/Semnet/stc/permutation/connectivity/parcellated/' # root


if not os.path.exists(out_path):
    os.makedirs(out_path)
subjects_dir = data_path 
print sys.argv
subject_inds = []#
for ss in sys.argv[1:]:
    subject_inds.append( int( ss ) )

subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]
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
            '/meg16_0123/160708/', #23 LD
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
            'MRI_meg16_0123' ,#23
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
wparc=1
parc_names=['aparc','aparc_mod','aparc2009','aparc2009_mod','RG']
labels_pre = mne.read_labels_from_annot('fsaverage', parc=parc_name[wparc], subjects_dir=subjects_dir)
path_ypos=results_path+ypos_names[wparc]
ltc=scipy.io.loadmat(path_ypos,mat_dtype=True); ypos_sort=np.squeeze(ltc['ypos_sort'])
labels=[labels_pre[jj] for jj in ypos_sort]
for lbl in labels:
    lbl.values.fill(1.0)
 
event_names=['Visual','Hear','Hand','Pwordc']
freq_names=['Theta','Alpha','Beta','Gamma']
#fname_vis='pacel_connmat_Visual_aparc_mod.mat'
#fname_hnd='pacel_connmat_Hand_aparc_mod.mat'

nsubs=len(list_all);
nnodes=int(parc_len[wparc])
grand_connmat=np.zeros((nnodes,nnodes,4,3,len(event_names),nsubs))#nnodes, nnodes, freq, time, nevent, nsubs
#vis_gamma_connmat=zeros(74,74,nsubs);
#hnd_alpha_connmat=zeros(74,74,nsubs);
#hnd_gamma_connmat=zeros(74,74,nsubs);

#% hnd_connmat=zeros(74,74,2,nsubs);
timesmin=['50ms','150ms','250ms']
timesmax=['250ms','350ms','450ms']
#thist=2
thisfreq=1#gamma
#thisev=1
refev=0
#thist=1
#N = 3; M = 3;                  %# grid size
#CONNECTED = 8;                 %# 4-/8- connected points
#
#%# which distance function
#if CONNECTED == 4,     distFunc = 'cityblock';
#elseif CONNECTED == 8, distFunc = 'chebychev'; end
#
#%# compute adjacency matrix
#nnodes=3
xv,yv = np.meshgrid(range(nnodes),range(nnodes));
xv=np.ravel(xv); yv=np.ravel(yv)# = X(:); Y = Y(:);
zv=np.vstack((yv,xv)).transpose()
adj = squareform( pdist(zv, metric='chebyshev') == 1 );
conn_sparse=sparse.csr_matrix(adj)
n_subjects=20
p_threshold=0.05
pthr=0.05
nlevels=4
#t_threshold = -stats.distributions.t.ppf(p_threshold/2., n_subjects - 1)
n_permutations=100
tail=0
freq1=1#alpha
freq2=3#gamma
factor_levels = [2,2]  # number of levels in each factor
effects = 'A:B'  # this is the default signature for computing all effects
return_pvals = True  
plt.close("all")

##    max_step=1;
t_threshold = mne.stats.f_threshold_mway_rm(n_subjects, factor_levels, effects, p_threshold)
for thist in range(3):
    for ii, meg in enumerate(ll):
#        print ii
        for evcnt, event in enumerate(event_names):
            fnameev_full=data_path+meg+'pacel_connmat_'+event+'_'+parc_names[wparc]+'.mat'
            ltc=scipy.io.loadmat(fnameev_full,mat_dtype=True); connmat=np.squeeze(ltc['connmat'])            
    
            for jj in range(connmat.shape[2]):#connmat nnodes, nnodes, freq, time
                for kk in range(connmat.shape[3]):
                    aa=np.squeeze(connmat[:,:,jj,kk])
                    aa=aa+aa.copy().transpose()+np.eye((aa.shape[0]))
                    connmat[:,:,jj,kk]=aa.copy()
            grand_connmat[:,:,:,:,evcnt,ii]=connmat.copy()
    MT=np.zeros((nnodes,nnodes,2))
    MTposb=np.zeros((nnodes,nnodes,2))
    for thisev in np.array([1,2]):     
        thisX=np.zeros((nnodes,nnodes,n_subjects,nlevels))#grand_connmat[:,:,thisfreq,thist,thisev,:]-grand_connmat[:,:,thisfreq,thist,refev,:]
        thisX[:,:,:,0]=grand_connmat[:,:,freq1,thist,thisev,:].copy()
        thisX[:,:,:,1]=grand_connmat[:,:,freq2,thist,thisev,:].copy()
        thisX[:,:,:,2]=grand_connmat[:,:,freq1,thist,refev,:].copy()
        thisX[:,:,:,3]=grand_connmat[:,:,freq2,thist,refev,:].copy()
    
        thisX=np.transpose(thisX,(2,3,0,1))#subxconxnodexnode
           
        for uvcnt in range(nnodes):
    #        print uvcnt
    #        MT[:,uvcnt]= mne.stats.ttest_1samp_no_p(thisX[:,:,uvcnt], sigma=1e-3, method='relative')
            MT[:,uvcnt,thisev-1], p= mne.stats.f_mway_rm(thisX[:,:,:,uvcnt], factor_levels=factor_levels,effects=effects, return_pvals=return_pvals,alpha=p_threshold)
    #    h, pc=mne.stats.fdr_correction(p)
        np.fill_diagonal(MT[:,:,thisev-1],0)
        #    MTabs=np.abs(MT)
        #    MTabs[MTabs<t_threshold]=0
        MTpos=MT[:,:,thisev-1].copy()
        if thisev==1:
            MTpos[np.abs(MTpos)<t_threshold]=0
            MTpos[np.abs(MTpos)>=t_threshold]=1
        else:
            MTpos[np.abs(MTpos)<t_threshold]=0
            MTpos[np.abs(MTpos)>=t_threshold]=-2
        MTposb[:,:,thisev-1]=MTpos.copy()
        
#        MTpos_1side=MT.copy()
#        MTpos_1side[MTpos_1side<t_threshold]=0
#        Tpos_sum=np.sum(MTpos_1side,1)
#        print Tpos_sum
#        out_fname=out_path+'MTpos_interaction_'+parc_names[wparc]+'_'+event_names[thisev]+'_'+event_names[refev]+'_'+freq_names[thisfreq]+'_'+timesmin[thist]+'_'+timesmax[thist]+'.mat'
#        scipy.io.savemat(out_fname,{'MTpos':MTpos})
#    comps_pos, comp_sizes_pos=get_components(MTpos)
#    MTpos=MTpos-np.eye(MTpos.shape[0])
#    cu_pos=np.unique(comps_pos)
#    Tpos_sum=cu_pos.copy()
#    for tci,thiscomp_pos in enumerate(cu_pos):
#        thisMT=MTpos_1side[comps_pos==thiscomp_pos,:]
#        Tpos_sum[tci]=np.sum(np.sum(thisMT,1))/2.
#    Tpos_max_perm=np.zeros((n_permutations,))    
#    for pcnt in range(n_permutations):
##        print pcnt
#        thisX=np.zeros((nnodes,nnodes,n_subjects,nlevels))#grand_connmat[:,:,thisfreq,thist,thisev,:]-grand_connmat[:,:,thisfreq,thist,refev,:]
#        rp1=np.random.permutation(n_subjects)
#        rp2=np.random.permutation(n_subjects)
#        thisc=rp1[:rp2[0]]
#        refc=np.setdiff1d(np.arange(n_subjects),thisc)
#        thisX[:,:,thisc,0]=grand_connmat[:,:,freq1,thist,thisev,thisc].copy()
#        thisX[:,:,refc,0]=grand_connmat[:,:,freq1,thist,refev,refc].copy()
#        thisX[:,:,thisc,1]=grand_connmat[:,:,freq2,thist,thisev,thisc].copy()
#        thisX[:,:,refc,1]=grand_connmat[:,:,freq2,thist,refev,refc].copy()
#        
#        thisX[:,:,thisc,2]=grand_connmat[:,:,freq1,thist,refev,thisc].copy()
#        thisX[:,:,refc,2]=grand_connmat[:,:,freq1,thist,thisev,refc].copy()
#        thisX[:,:,thisc,3]=grand_connmat[:,:,freq2,thist,refev,thisc].copy()
#        thisX[:,:,refc,3]=grand_connmat[:,:,freq2,thist,thisev,refc].copy()
#        
##        thisX[:,:,:,0]=grand_connmat[:,:,freq1,thist,thisev,:].copy()
##        thisX[:,:,:,1]=grand_connmat[:,:,freq2,thist,thisev,:].copy()
##        thisX[:,:,:,2]=grand_connmat[:,:,freq1,thist,refev,:].copy()
##        thisX[:,:,:,3]=grand_connmat[:,:,freq2,thist,refev,:].copy()
#        
#    
#        thisX=np.transpose(thisX,(2,3,0,1))#subxconxnodexnode
#        MT=np.zeros((nnodes,nnodes))    
#        for uvcnt in range(nnodes):
#    #        print uvcnt
#    #        MT[:,uvcnt]= mne.stats.ttest_1samp_no_p(thisX[:,:,uvcnt], sigma=1e-3, method='relative')
#            MT[:,uvcnt], p= mne.stats.f_mway_rm(thisX[:,:,:,uvcnt], factor_levels=factor_levels,effects=effects, return_pvals=return_pvals,alpha=p_threshold)
#        np.fill_diagonal(MT,0)
#    #    h, pc=mne.stats.fdr_correction(p)
#   
##        thisX=grand_connmat[:,:,thisfreq,thist,thisev,:]-grand_connmat[:,:,thisfreq,thist,refev,:]
##        thisX=np.transpose(thisX,(2,0,1))
##        MT=np.zeros((nnodes,nnodes))  
##        signvect=np.sign(np.random.randn(n_subjects))
##        signmat=np.tile(signvect,(nnodes,1)).transpose()
##        for uvcnt in range(thisX.shape[2]):
##    #            print uvcnt
##            MT[:,uvcnt]= mne.stats.ttest_1samp_no_p(signmat*thisX[:,:,uvcnt], sigma=1e-3, method='relative')
##        MTabs=np.abs(MT)
##        MTabs[MTabs<t_threshold]=0
#        MTpos_perm=MT.copy()
#        MTpos_perm[MTpos_perm<t_threshold]=0
#        MTpos_perm[MTpos_perm>=t_threshold]=1
#        #MTpos=MTpos-np.eye(MTpos.shape[0])
#        MTpos_1side=MT.copy()
#        MTpos_1side[MTpos_1side<t_threshold]=0
#        Tpos_sum_perm=np.sum(MTpos_1side,1)
##        comps_pos_perm, comp_sizes_pos=get_components(MTpos_perm)
##        MTpos_perm=MTpos_perm-np.eye(MTpos_perm.shape[0])
##        cu_pos_perm=np.unique(comps_pos_perm)
##        Tpos_sum_perm=cu_pos_perm.copy()
##        for tci,thiscomp_pos in enumerate(cu_pos_perm):
##            thisMT=MTpos_1side[comps_pos_perm==thiscomp_pos,:]
##            Tpos_sum_perm[tci]=np.sum(np.sum(thisMT,1))/2.
#        Tpos_max_perm[pcnt]=np.max(Tpos_sum_perm)
#    pvals=np.zeros((len(Tpos_sum),))
#    for pvcnt in range(len(pvals)):
#        pvals[pvcnt]=np.sum(Tpos_max_perm>=Tpos_sum[pvcnt])/np.float(n_permutations)
#    print pvals
#    pval_sig=np.where(pvals<pthr)[0]
#    this_MTpos=np.zeros(MTpos.shape)
#    if len(pval_sig)>0:
##        for pvs,pvalsig in enumerate(pval_sig):
##        thisp=pvals[pvalsig]
#        this_MTpos[pval_sig,:]=MTpos[pval_sig,:].copy()    
#        out_fname=out_path+'NBSnetwork_interaction_pos_'+parc_names[wparc]+'_'+event_names[thisev]+'_'+event_names[refev]+'_'+freq_names[thisfreq]+'_'+timesmin[thist]+'_'+timesmax[thist]+'_'+str(p_threshold)[2:]+'.mat'
#        scipy.io.savemat(out_fname,{'this_MTpos':this_MTpos})
    
#            
    
    ###NEG!
    
    
#    MTneg=MT.copy()
#    MTneg[MTneg>-t_threshold]=0
#    MTneg[MTneg<=-t_threshold]=1
#    #MTpos=MTpos-np.eye(MTpos.shape[0])
#    MTneg_1side=MT.copy()
#    MTneg_1side[MTneg_1side>-t_threshold]=0
#    comps_neg, comp_sizes_neg=get_components(MTneg)
#    MTneg=MTneg-np.eye(MTneg.shape[0])
#    cu_neg=np.unique(comps_neg)
#    Tneg_sum=cu_neg.copy()
#    for tci,thiscomp_neg in enumerate(cu_neg):
#        thisMT=MTneg_1side[comps_neg==thiscomp_neg,:]
#        Tneg_sum[tci]=np.sum(np.sum(thisMT,1))/2.
#    Tneg_max_perm=np.zeros((n_permutations,))    
#    for pcnt in range(n_permutations):
##        print pcnt        
#        thisX=grand_connmat[:,:,thisfreq,thist,thisev,:]-grand_connmat[:,:,thisfreq,thist,refev,:]
#        thisX=np.transpose(thisX,(2,0,1))
#        MT=np.zeros((nnodes,nnodes))  
#        signvect=np.sign(np.random.randn(n_subjects))
#        signmat=np.tile(signvect,(nnodes,1)).transpose()
#        for uvcnt in range(thisX.shape[2]):
#    #            print uvcnt
#            MT[:,uvcnt]= mne.stats.ttest_1samp_no_p(signmat*thisX[:,:,uvcnt], sigma=1e-3, method='relative')
#        MTabs=np.abs(MT)
#        MTabs[MTabs<t_threshold]=0
#        MTneg_perm=MT.copy()
#        MTneg_perm[MTneg_perm>-t_threshold]=0
#        MTneg_perm[MTneg_perm<=-t_threshold]=1
#        #MTneg=MTneg-np.eye(MTneg.shape[0])
#        MTneg_1side=MT.copy()
#        MTneg_1side[MTneg_1side>-t_threshold]=0
#        comps_neg_perm, comp_sizes_neg=get_components(MTneg_perm)
#        MTneg_perm=MTneg_perm-np.eye(MTneg_perm.shape[0])
#        cu_neg_perm=np.unique(comps_neg_perm)
#        Tneg_sum_perm=cu_neg_perm.copy()
#        for tci,thiscomp_neg in enumerate(cu_neg_perm):
#            thisMT=MTneg_1side[comps_neg_perm==thiscomp_neg,:]
#            Tneg_sum_perm[tci]=np.sum(np.sum(thisMT,1))/2.
#        Tneg_max_perm[pcnt]=np.min(Tneg_sum_perm)
#    pvals=np.zeros((len(Tneg_sum),))
#    for pvcnt in range(len(pvals)):
#        pvals[pvcnt]=np.sum(np.abs(Tneg_max_perm)>=np.abs(Tneg_sum[pvcnt]))/np.float(n_permutations)
#    print pvals
#    pval_sig=np.where(pvals<pthr)[0]
#    if len(pval_sig)>0:
#        for pvs,pvalsig in enumerate(pval_sig):
#            this_MTneg=MTneg.copy()
#            this_MTneg[np.where(comps_neg!=cu_neg[pvalsig])[0],:]=0
#            out_fname=out_path+'NBSnetwork_neg_cluster'+str(pvs)+'_'+parc_names[wparc]+'_'+event_names[thisev]+'_'+event_names[refev]+'_'+freq_names[thisfreq]+'_'+timesmin[thist]+'_'+timesmax[thist]+'_'+str(p_threshold)[2:]+'.mat'
#            scipy.io.savemat(out_fname,{'this_MTneg':this_MTneg})
#    
    MTposf=MTposb[:,:,0]+MTposb[:,:,1]
    MTposf_bin=MTposf.copy()
    MTposf_bin[MTposf_bin!=0]=1
    ss=np.sum(MTposf_bin,1)
    print ss
    sshub=np.where(ss>np.mean(ss)+2*np.std(ss))[0]
    ssnohub=np.setdiff1d(np.arange(len(labels)),sshub)
    MTposf[ssnohub,:]=0
    fig = plt.figure(num=thist, figsize=(8, 6), facecolor='black')
    label_names = [label.name for label in labels]
    label_colorsp = [label.color for label in labels]
    label_colors=[(lc[0],lc[1],lc[2],1.0) for lc in label_colorsp]
    
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
    #            plt.figure()#close("all")
    #            con_res_mean=np.squeeze(np.mean(thispval,axis=2))
    
    plot_connectivity_circle(MTposf, range(len(label_names)),colormap='jet',n_lines=None,
                         node_angles=node_angles, node_colors=label_colors,colorbar=True,linewidth=0.5,#, fontsize_title=14.,fontsize_names=6.,fontsize_colorbar=12.,colorbar_size=0.4,
                         fig=fig, padding=1.0)#subplot=(1, 5, which_parc+1),
#                                 
                                 
                                 
#MTneg=MT.copy()
#MTneg[MTneg>-t_threshold]=0
#MTneg[MTneg<=-t_threshold]=1
#comps_neg, comp_sizes_neg=get_components(MTneg)


    
#    t,c,cp,h=mne.stats.permutation_cluster_1samp_test(thisX, threshold=t_threshold, n_permutations=n_permutations, 
#    tail=tail,  t_power=1,n_jobs=4,connectivity=conn_sparse)#, step_down_p=0.05,max_step=2,
#    for n1 in range(nnodes):
#        for n2 in range(n1,nnodes):
#        grand_connmat[:,:,:,:,evcnt,ii]=connmat.copy()
#for evcnt
#        vis_gamma_connmat(:,:,ii)=squeeze(connmat(:,:,4,thist));
#        
#        fnamehnd_full=[data_path,list_all{ii},fname_hnd];
#        load(fnamehnd_full)
#        for jj=1:size(connmat,3)%freqs
#            for kk=1:size(connmat,4)%times
#                aa=squeeze(connmat(:,:,jj,kk));
#                aa=aa+aa'+eye(size(aa));
#                connmat(:,:,jj,kk)=aa;
#            end
#        end
#        hnd_alpha_connmat(:,:,ii)=squeeze(connmat(:,:,2,thist));
#        hnd_gamma_connmat(:,:,ii)=squeeze(connmat(:,:,4,thist));
#
#grand_connmat=cat(3,vis_alpha_connmat,vis_gamma_connmat,hnd_alpha_connmat,hnd_gamma_connmat);
#fname_out=[data_path,'/NBS1.2/Rezvan/','Matrix_vishnd_alphagamma_',timesmin{thist},'_',timesmax{thist},'.mat'];
#save(fname_out,'grand_connmat')
#gconnmat2=cat(3,vis_gamma_connmat,hnd_gamma_connmat);
#fname_out=[data_path,'/NBS1.2/Rezvan/','Matrix_vishnd_gamma_',timesmin{thist},'_',timesmax{thist},'.mat'];
#save(fname_out,'gconnmat2')
#
#gconnmat3=cat(3,vis_alpha_connmat,hnd_alpha_connmat);%vis_alpha_connmat-hnd_alpha_connmat;%
#fname_out=[data_path,'/NBS1.2/Rezvan/','Matrix_vishnd_alpha_',timesmin{thist},'_',timesmax{thist},'.mat'];
#save(fname_out,'gconnmat3')
#
#
#dmat=cat(1,eye(nsubs),eye(nsubs));
#dmat=[dmat,[ones(nsubs,1);-1*ones(nsubs,1)]];
#% dmat=ones(nsubs,1);
#fname_out=[data_path,'/NBS1.2/Rezvan/','designMatrix_vishnd_gamma_',timesmin{thist},'_',timesmax{thist},'.mat'];
#save(fname_out,'dmat')
#load('/imaging/rf02/Semnet/NBS1.2/Rezvan/nodeinfo_parc1.mat')
#nodenames={};
#for nc=1:size(centre_pos,1)
#    if mod(nc,2)==1 %left
#    bb=['L_',num2str(((nc+1)/2))];
#    else
#    bb=['R_',num2str((nc/2))];
#    end
#    nodenames{nc,1}=bb;
#end
#
#cont=vis_alpha_connmat-hnd_alpha_connmat;
#nnodes=size(centre_pos,1);
#pval=zeros(nnodes,nnodes);
#for p1=1:nnodes
#    for p2=1:nnodes
#       [h, pval(p1,p2)]=ttest(squeeze(cont(p1,p2,:)));
#    end
#end
#
#
#    
    
    
    

    
    
    
    
    
    
    
    
    
    