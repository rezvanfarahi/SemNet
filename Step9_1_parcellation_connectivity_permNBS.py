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
 
event_names=['Visual','Hear','Hand','Pwordc']
#fname_vis='pacel_connmat_Visual_aparc_mod.mat'
#fname_hnd='pacel_connmat_Hand_aparc_mod.mat'

nsubs=len(list_all);
nnodes=74
grand_connmat=np.zeros((nnodes,nnodes,4,3,len(event_names),nsubs))#nnodes, nnodes, freq, time, nevent, nsubs
#vis_gamma_connmat=zeros(74,74,nsubs);
#hnd_alpha_connmat=zeros(74,74,nsubs);
#hnd_gamma_connmat=zeros(74,74,nsubs);

#% hnd_connmat=zeros(74,74,2,nsubs);
timesmin=['50ms','150ms','250ms']
timesmax=['250ms','350ms','450ms']
#thist=2
thisfreq=1#gamma
thisev=0
refev=3
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
t_threshold = -stats.distributions.t.ppf(p_threshold/2., n_subjects - 1)
n_permutations=10000
tail=0
for ii, meg in enumerate(ll):
    print ii
    for evcnt, event in enumerate(event_names):
        fnameev_full=data_path+meg+'pacel_connmat_'+event+'_aparc_mod.mat'
        ltc=scipy.io.loadmat(fnameev_full,mat_dtype=True); connmat=np.squeeze(ltc['connmat'])            

        for jj in range(connmat.shape[2]):#connmat nnodes, nnodes, freq, time
            for kk in range(connmat.shape[3]):
                aa=np.squeeze(connmat[:,:,jj,kk])
                aa=aa+aa.copy().transpose()+np.eye((aa.shape[0]))
                connmat[:,:,jj,kk]=aa.copy()
        grand_connmat[:,:,:,:,evcnt,ii]=connmat.copy()
for thist in np.array([1]):
    thisX=grand_connmat[:,:,thisfreq,thist,thisev,:]-grand_connmat[:,:,thisfreq,thist,refev,:]
    thisX=np.transpose(thisX,(2,0,1))
    
    t,c,cp,h=mne.stats.permutation_cluster_1samp_test(thisX, threshold=t_threshold, n_permutations=n_permutations, 
    tail=tail,  t_power=1,n_jobs=4,connectivity=conn_sparse)#, step_down_p=0.05,max_step=2,
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
    
    
    

    
    
    
    
    
    
    
    
    
    