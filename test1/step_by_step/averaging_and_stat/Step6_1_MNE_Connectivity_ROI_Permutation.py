# -*- coding: utf-8 -*-
"""
Created on Sat Oct 29 01:28:31 2016

@author: rf02
"""
print(__doc__)
import sys

#sys.path.append('/imaging/local/software/python_packages/nibabel/1.3.0')
#sys.path.append('/imaging/local/software/python_packages/pysurfer/v0.4')
#sys.path.append('/imaging/local/software/python_packages/scikit-learn/v0.14.1')
#sys.path.append('/imaging/local/software/mne_python/latest')

# for qsub
# (add ! here if needed)
sys.path.append('/imaging/local/software/EPD/latest/x86_64/bin/python')
#sys.path.insert(1,'/imaging/rf02/TypLexMEG/mne_python0.9')
sys.path.insert(1,'/imaging/rf02/mne_python_11')

#sys.path.insert(1,'/imaging/local/software/mne_python/latest_v0.9')
#sys.path.append('/imaging/local/software/python_packages/scikit-learn/v0.16.1')
sys.path.insert(1,'/imaging/local/software/python_packages/scikit-learn/v0.13.1/lib.linux-x86_64-2.7/sklearn/externals')
sys.path.append('/home/rf02/rezvan/test1')
#sys.path.append('/imaging/local/software/python_packages/scikit-learn/v0.14.1')
#sys.path.insert(1,'/imaging/local/software/python_packages/scikit-learn/v0.16.1')
sys.path.insert(1,'/imaging/rf02/scikit-learn-0.16.1')
sys.path.insert(1,'/home/rf02/.local/lib/python2.7/site-packages')

import joblib
###

import os
import numpy as np
import mne
import scipy
#import sklearn

from mne.stats import permutation_t_test, ttest_1samp_no_p, fdr_correction
from numpy import *
import numpy

in_path='/imaging/rf02/TypLexMEG/icaanalysis_results/ROI/concrete_coh_ppc_1stpc_bands_1stpc.mat';
SC_pca_allcb_pre=scipy.io.loadmat(in_path,mat_dtype=True); SC_pca_allcb=np.squeeze(SC_pca_allcb_pre['SC_pca_allcb'])

in_path='/imaging/rf02/TypLexMEG/icaanalysis_results/ROI/abstract_coh_ppc_1stpc_bands_1stpc.mat'; 
SC_pca_allab_pre=scipy.io.loadmat(in_path,mat_dtype=True); SC_pca_allab=np.squeeze(SC_pca_allab_pre['SC_pca_allab'])
SC_pca_allsb=np.subtract(SC_pca_allcb,SC_pca_allab)
aa=SC_pca_allsb[:,:,:,0];bb=aa[np.abs(aa)>0];
nsamp=SC_pca_allsb.shape[-1]
ntest=len(bb)
X=np.zeros((nsamp,ntest))
for ii in range(SC_pca_allsb.shape[-1]):
    aa=SC_pca_allsb[:,:,:,ii];bb=aa[np.abs(aa)>0];
    X[ii,:]=bb
    
t,p,h=permutation_t_test(X)    
t1=ttest_1samp_no_p(X)   
pval=scipy.stats.t.sf(np.abs(t1),16)*2
r,pc=fdr_correction(pval)#,alpha=0.5)# method='negcorr'
    
#in_path='/imaging/rf02/TypLexMEG/icaanalysis_results/ROI/concrete_coh_ppc_1stpc_bands_all.mat';
#in_path='/imaging/rf02/TypLexMEG/icaanalysis_results/ROI/abstract_coh_ppc_1stpc_bands_all.mat';
this_band=3
this_time=1
in_path='/imaging/rf02/TypLexMEG/icaanalysis_results/ROI/concrete_coh_ppc_1stpc_bands_all.mat';
SC_pca_allcb_pre=scipy.io.loadmat(in_path,mat_dtype=True); SC_pca_allc=np.squeeze(SC_pca_allcb_pre['SC_pca_allc'])
SC_pca_allc_tb=np.squeeze(SC_pca_allc[:,:,:,this_band,:])
in_path='/imaging/rf02/TypLexMEG/icaanalysis_results/ROI/abstract_coh_ppc_1stpc_bands_all.mat'; 
SC_pca_allab_pre=scipy.io.loadmat(in_path,mat_dtype=True); SC_pca_alla=np.squeeze(SC_pca_allab_pre['SC_pca_alla'])
SC_pca_alla_tb=np.squeeze(SC_pca_alla[:,:,:,this_band,:])

SC_pca_allsb=np.subtract(SC_pca_allc_tb,SC_pca_alla_tb)
aa=SC_pca_allsb[:,:,this_time,0];bb=aa[np.abs(aa)>0];
nsamp=SC_pca_allsb.shape[-1]
ntest=len(bb)
X=np.zeros((nsamp,ntest))
for ii in range(SC_pca_allsb.shape[-1]):
    aa=SC_pca_allsb[:,:,this_time,ii];bb=aa[np.abs(aa)>0];
    X[ii,:]=bb
    
t,p,h=permutation_t_test(X)    
t1=ttest_1samp_no_p(X)   
pval=scipy.stats.t.sf(np.abs(t1),16)*2
r,pc=fdr_correction(pval)#,alpha=0.5)# method='negcorr'
print p.min()
print pc.min()
#
in_path='/imaging/rf02/TypLexMEG/icaanalysis_results/ROI/concrete_mi_bands_1stpc.mat';
SC_pca_allcb_pre=scipy.io.loadmat(in_path,mat_dtype=True); SC_pca_allcb=np.squeeze(SC_pca_allcb_pre['SC_pca_allcb'])

in_path='/imaging/rf02/TypLexMEG/icaanalysis_results/ROI/abstract_mi_bands_1stpc.mat';
SC_pca_allab_pre=scipy.io.loadmat(in_path,mat_dtype=True); SC_pca_allab=np.squeeze(SC_pca_allab_pre['SC_pca_allab'])
SC_pca_allsb=np.subtract(SC_pca_allcb,SC_pca_allab)
aa=SC_pca_allsb[:,:,:,0];bb=aa[np.abs(aa)>0];
nsamp=SC_pca_allsb.shape[-1]
ntest=len(bb)
X=np.zeros((nsamp,ntest))
for ii in range(SC_pca_allsb.shape[-1]):
    aa=SC_pca_allsb[:,:,:,ii];bb=aa[np.abs(aa)>0];
    X[ii,:]=bb
    
t,p,h=permutation_t_test(X)    
t1=ttest_1samp_no_p(X)   
pval=scipy.stats.t.sf(np.abs(t1),16)*2
r,pc=fdr_correction(pval)

this_band=2
this_time=1
in_path='/imaging/rf02/TypLexMEG/icaanalysis_results/ROI/concrete_mi_bands_all.mat';
SC_pca_allc_pre=scipy.io.loadmat(in_path,mat_dtype=True); SC_pca_allc=np.squeeze(SC_pca_allc_pre['SC_pca_allc'])
SC_pca_allc_tb=np.squeeze(SC_pca_allc[:,:,:,this_band,:])
in_path='/imaging/rf02/TypLexMEG/icaanalysis_results/ROI/abstract_mi_bands_all.mat';
SC_pca_alla_pre=scipy.io.loadmat(in_path,mat_dtype=True); SC_pca_alla=np.squeeze(SC_pca_alla_pre['SC_pca_alla'])
SC_pca_alla_tb=np.squeeze(SC_pca_alla[:,:,:,this_band,:])

SC_pca_allsb=np.subtract(SC_pca_allc_tb,SC_pca_alla_tb)
aa=SC_pca_allsb[:,:,this_time,0];bb=aa[np.abs(aa)>0];
nsamp=SC_pca_allsb.shape[-1]
ntest=len(bb)
X=np.zeros((nsamp,ntest))
for ii in range(SC_pca_allsb.shape[-1]):
    aa=SC_pca_allsb[:,:,this_time,ii];bb=aa[np.abs(aa)>0];
    X[ii,:]=bb
    
t,p,h=permutation_t_test(X)
print p.min()    
t1=ttest_1samp_no_p(X)   
pval=scipy.stats.t.sf(np.abs(t1),16)*2
r,pc=fdr_correction(pval)
print pc.min()
#plt.hist(p,bins=50)
ss=np.mean(SC_pca_allsb,axis=-1)
plt.imshow(ss[:,:,0], cmap='coolwarm', aspect='equal', interpolation="nearest")