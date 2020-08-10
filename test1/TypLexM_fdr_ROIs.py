# -*- coding: utf-8 -*-
"""
Created on Mon Jun  1 23:29:48 2015

@author: rf02
"""

"""
=========================================================================
Compute source space connectivity and visualize it using a circular graph
=========================================================================

This example computes the all-to-all connectivity between 68 regions in
source space based on dSPM inverse solutions and a FreeSurfer cortical
parcellation. The connectivity is visualized using a circular graph which
is ordered based on the locations of the regions.
"""

# Authors: Martin Luessi <mluessi@nmr.mgh.harvard.edu>
#          Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
#          Nicolas P. Rougier (graph code borrowed from his matplotlib gallery)
#
# License: BSD (3-clause)

print(__doc__)
import sys

#sys.path.append('/imaging/local/software/python_packages/scikit-learn/v0.14.1')
sys.path.insert(1,'/imaging/local/software/mne_python/latest_v0.9')
sys.path.insert(1,'/imaging/local/software/mne_python/latest_v0.9/scikit-learn')

# for qsub
# (add ! here if needed) 
sys.path.append('/imaging/local/software/anaconda/latest/x86_64/bin/python')
#sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
###

import numpy as np
import mne
from mne.io import Raw
from mne.minimum_norm import (read_inverse_operator, apply_inverse, apply_inverse_epochs)
import scipy.io
from mne import find_events
from mne.connectivity import (spectral_connectivity, phase_slope_index)
import numpy as np
from mne import read_labels_from_annot
import matplotlib.pyplot as plt


data_path = '/imaging/rf02/TypLexMEG/' # root directory for your MEG data
subjects_dir = '/imaging/rf02/TypLexMEG/'    # where your MRI subdirectories are
out_path = '/imaging/rf02/TypLexMEG/icaanalysis_results/stc/UVttest/time_resolved/' # root directory for your MEG data
event_path = event_path = '/group/erp/data/olaf.hauk/MEG/TypLexM/data/MF22_new/'   # where event files are

label_path = '/imaging/rf02/TypLexMEG/fsaverage/label/'

inv_path = '/imaging/rf02/TypLexMEG/'
inv_fname = 'InverseOperator_fft_1_48_clean_ica_EMEG-inv.fif'

# get indices for subjects to be processed from command line input
# 
print sys.argv
subject_inds = []
for ss in sys.argv[1:]:
   subject_inds.append( int( ss ) )



subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]#0,1,2,3,4,5,6,7,8]#,
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
trange=np.arange(5,550,5)

print "ll:"
print ll
Dmat_cncrt=np.zeros((17,64,4,1201))
Dmat_abs=np.zeros((17,64,4,1201))
Dmat_sub=np.zeros((17,64,4,1201))
Dmat_cncrt_r=np.zeros((17,64,4,len(trange)))
Dmat_abs_r=np.zeros((17,64,4,len(trange)))
Dmat_sub_r=np.zeros((17,64,4,len(trange)))
for cnt, meg in enumerate(ll):
    
    print cnt
    fname_cncrt = data_path + meg + 'firstMorphed_SemLoc_icaclean_Concrete_Source_Evoked_m500_700'
    stc_cncrt = mne.read_source_estimate(fname_cncrt)
    path_c=data_path+ meg+'degree_proportional_25_pca_ppc_concrete.mat'
    ltc=scipy.io.loadmat(path_c,mat_dtype=True); deg_c=ltc['deg_c']
    Dmat_cncrt[cnt,:,:,:]=deg_c
    path_a=data_path+ meg+'degree_proportional_25_pca_ppc_abstract.mat'
    lta=scipy.io.loadmat(path_a,mat_dtype=True); deg_a=lta['deg_a']
    Dmat_abs[cnt,:,:,:]=deg_a
    
    deg_s=deg_c-deg_a;
    Dmat_sub[cnt,:,:,:]=deg_s
for cc,bc in enumerate(trange):
    b1=bc; b2=b1+10-1
    Dmat_cncrt_r[:,:,:,cc]=np.subtract(Dmat_cncrt[:,:,:,b1:b2].mean(axis=3),Dmat_cncrt[:,:,:,200:300].mean(axis=3))
    Dmat_abs_r[:,:,:,cc]=np.subtract(Dmat_abs[:,:,:,b1:b2].mean(axis=3),Dmat_abs[:,:,:,200:300].mean(axis=3))
    Dmat_sub_r[:,:,:,cc]=np.subtract(Dmat_cncrt[:,:,:,b1:b2] .mean(axis=3),Dmat_abs[:,:,:,b1:b2] .mean(axis=3))
print "resample fin"    
MTc=np.zeros(Dmat_cncrt_r.shape[1:])
MTa=np.zeros(Dmat_abs_r.shape[1:])
MTs=np.zeros(Dmat_sub_r.shape[1:])
MPc=np.zeros(Dmat_cncrt_r.shape[1:])
MPa=np.zeros(Dmat_abs_r.shape[1:])
MPs=np.zeros(Dmat_sub_r.shape[1:])
for ii in range(MTc.shape[0]):
    print ii
    for jj in range(MTc.shape[1]):
        for kk in range(MTc.shape[2]):
            MTc[ii,jj,kk],MPc[ii,jj,kk],h=mne.stats.permutation_t_test(Dmat_cncrt_r[:,ii:ii+1,jj,kk], tail=0,n_permutations=1024)
            MTa[ii,jj,kk],MPa[ii,jj,kk],h=mne.stats.permutation_t_test(Dmat_abs_r[:,ii:ii+1,jj,kk], tail=0,n_permutations=1024)
            MTs[ii,jj,kk],MPs[ii,jj,kk],h=mne.stats.permutation_t_test(Dmat_sub_r[:,ii:ii+1,jj,kk], tail=0,n_permutations=1024)
print "permutation done"

#for jj in range(MTc.shape[1]):
#    for kk in range(MTc.shape[2]):
#        MTc[:,jj,kk],MPc[:,jj,kk],h=mne.stats.permutation_t_test(Dmat_cncrt_r[:,:,jj,kk], tail=0,n_permutations=1024)
#        MTa[:,jj,kk],MPa[:,jj,kk],h=mne.stats.permutation_t_test(Dmat_abs_r[:,:,jj,kk], tail=0,n_permutations=1024)
#        MTs[:,jj,kk],MPs[:,jj,kk],h=mne.stats.permutation_t_test(Dmat_sub_r[:,:,jj,kk], tail=0,n_permutations=1024)
#print "permutation done"


#MTc,MPc,h=mne.stats.permutation_t_test(Dmat_cncrt_r.reshape((Dmat_cncrt_r.shape[0],Dmat_cncrt_r.shape[1]*Dmat_cncrt_r.shape[2]*Dmat_cncrt_r.shape[3])), tail=0,n_permutations=1024)
#MTc=MTc.reshape((Dmat_cncrt_r.shape[1],Dmat_cncrt_r.shape[2],Dmat_cncrt_r.shape[3]))
#MPc=MPc.reshape((Dmat_cncrt_r.shape[1],Dmat_cncrt_r.shape[2],Dmat_cncrt_r.shape[3]))
#
#MTa,MPa,h=mne.stats.permutation_t_test(Dmat_abs_r.reshape((Dmat_abs_r.shape[0],Dmat_abs_r.shape[1]*Dmat_abs_r.shape[2]*Dmat_abs_r.shape[3])), tail=0,n_permutations=1024)
#MTa=MTa.reshape((Dmat_abs_r.shape[1],Dmat_abs_r.shape[2],Dmat_abs_r.shape[3]))
#MPa=MPa.reshape((Dmat_abs_r.shape[1],Dmat_abs_r.shape[2],Dmat_abs_r.shape[3]))
#
#MTs,MPs,h=mne.stats.permutation_t_test(Dmat_sub_r.reshape((Dmat_sub_r.shape[0],Dmat_sub_r.shape[1]*Dmat_sub_r.shape[2]*Dmat_sub_r.shape[3])), tail=0,n_permutations=1024)
#MTs=MTs.reshape((Dmat_sub_r.shape[1],Dmat_sub_r.shape[2],Dmat_sub_r.shape[3]))
#MPs=MPs.reshape((Dmat_sub_r.shape[1],Dmat_sub_r.shape[2],Dmat_sub_r.shape[3]))

print "permutation done"


####FDR correction
#pvalc=MPc#np.zeros(MTc.shape)
#pvala=MPa#np.zeros(MTa.shape)
#pvals=MPs#.zeros(MTs.shape)
#pvalc_c=np.zeros(MTc.shape)
#pvala_c=np.zeros(MTa.shape)
#pvals_c=np.zeros(MTs.shape)
##
#for ii in range(pvalc.shape[1]):
#    for jj in range(pvalc.shape[2]):
#        #pvalc[:,ii,jj]=scipy.stats.t.sf(np.abs(MTc[:,ii,jj]),16)*2
#        hh,pvalc_c[:,ii,jj]=mne.stats.fdr_correction(pvalc[:,ii,jj])
#        #pvala[:,ii,jj]=scipy.stats.t.sf(np.abs(MTa[:,ii,jj]),16)*2
#        hh,pvala_c[:,ii,jj]=mne.stats.fdr_correction(pvala[:,ii,jj])
#        #pvals[:,ii,jj]=scipy.stats.t.sf(np.abs(MTs[:,ii,jj]),16)*2
#        hh,pvals_c[:,ii,jj]=mne.stats.fdr_correction(pvals[:,ii,jj])
#print pvals_c[pvals_c<0.05]
#print pvalc_c[pvalc_c<0.05]
#print pvala_c[pvala_c<0.05]

pvalc_c=MPc.copy()
Pvalc=-10*np.log10(pvalc_c)
Pvalc[Pvalc<-10*np.log10(0.05)]=0
Pvalc[MTc<0]=-Pvalc[MTc<0]

pvala_c=MPa.copy()
Pvala=-10*np.log10(pvala_c)
Pvala[Pvala<-10*np.log10(0.05)]=0
Pvala[MTa<0]=-Pvala[MTa<0]

pvals_c=MPs.copy()
Pvals=-10*np.log10(pvals_c)
Pvals[Pvals<-10*np.log10(0.05)]=0
Pvals[MTs<0]=-Pvals[MTs<0]


      
fsave_vertices = [np.arange(10242), np.arange(10242)]
Xc=np.zeros((20484, Pvalc.shape[2]))
Xa=np.zeros((20484, Pvala.shape[2]))
Xs=np.zeros((20484, Pvals.shape[2]))


print "deg files loaded"
labelss = mne.read_labels_from_annot(subject='fsaverage', parc='rezvan_annot', subjects_dir=data_path)
labelss.remove(labelss[-2])
labelss.remove(labelss[-1])
for bn,band in zip(range(Pvalc.shape[1]),['Theta','Alpha','Beta','Gamma']):
    print band
    for ccc,label in enumerate(labelss):  
        print ccc
        bb=stc_cncrt.in_label(label)
        if ccc%2==0:
            for cnl in bb.lh_vertno:
                Xc[cnl,:]=Pvalc[ccc,bn,:]
                Xa[cnl,:]=Pvala[ccc,bn,:]
                Xs[cnl,:]=Pvals[ccc,bn,:]
    
        else: 
            for cnr in bb.rh_vertno:
                Xc[cnr+10242,:]=Pvalc[ccc,bn,:]
                Xa[cnr+10242,:]=Pvala[ccc,bn,:]
                Xs[cnr+10242,:]=Pvals[ccc,bn,:]
                
        
       
    tmin1=5
    tstep1=5
    vertices_to = [np.arange(10242), np.arange(10242)]
    Xc_stc = mne.SourceEstimate(Xc, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
    Xa_stc = mne.SourceEstimate(Xa, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
    Xs_stc = mne.SourceEstimate(Xs, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
    
    out_file_Xc=out_path + 'Perttest_10logPval_firstMorphed_SemLoc_icaclean_degree_proportional_25_pca_ppc_concretebase_meanresampled_5_550_10ms_5overlap_'+band
    Xc_stc.save(out_file_Xc)
    
    out_file_Xa=out_path + 'Perttest_10logPval_firstMorphed_SemLoc_icaclean_degree_proportional_25_pca_ppc_abstractbase_meanresampled_5_550_10ms_5overlap_'+band
    Xa_stc.save(out_file_Xa)
    
    out_file_Xs=out_path  + 'Perttest_10logPval_firstMorphed_SemLoc_icaclean_degree_proportional_25_pca_ppc_subtract_meanresampled_5_550_10ms_5overlap_'+band
    Xs_stc.save(out_file_Xs)
#    

