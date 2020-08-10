# -*- coding: utf-8 -*-
"""
Created on Sun Jun 28 03:16:38 2015

@author: rf02
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 13:42:14 2015

@author: rf02
"""

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
trange=np.arange(1,10,2)#np.arange(10,91,20)
graphpath='icaanalysis_results/jpg/time_resolved/ROI_connectivity/graph/'

print "ll:"
print ll
Dmat_cncrt=np.zeros((17,46,4,2))
Dmat_abs=np.zeros((17,46,4,2))
Dmat_sub=np.zeros((17,46,4,2))
Dmat_cncrt_r=np.zeros((17,46,4,len(trange)))
Dmat_abs_r=np.zeros((17,46,4,len(trange)))
Dmat_sub_r=np.zeros((17,46,4,len(trange)))
for cnt, meg in enumerate(ll):
    
    print cnt
    fname_cncrt = data_path + meg + 'firstMorphed_SemLoc_icaclean_Concrete_Source_Evoked_m500_700'
    stc_cncrt = mne.read_source_estimate(fname_cncrt)
    path_c=data_path+ meg+'degrees_proportional_15_pca_coh_concrete_resampled_2win'
    ltc=scipy.io.loadmat(path_c,mat_dtype=True); deg_c=ltc['deg_cr']
    degc=deg_c.copy()
    for dii in range(deg_c.shape[1]):
        for djj in range(deg_c.shape[2]):
            degc[:,dii,djj]=(deg_c[:,dii,djj]-deg_c[:,dii,djj].mean())/np.std(deg_c[:,dii,djj])
    Dmat_cncrt[cnt,:,:,:]=deg_c
    path_a=data_path+ meg+'degrees_proportional_15_pca_coh_abstract_resampled_2win'
    lta=scipy.io.loadmat(path_a,mat_dtype=True); deg_a=lta['deg_ar']
    dega=deg_a.copy()
    for dii in range(deg_a.shape[1]):
        for djj in range(deg_a.shape[2]):
            dega[:,dii,djj]=(deg_a[:,dii,djj]-deg_a[:,dii,djj].mean())/np.std(deg_a[:,dii,djj])
    Dmat_abs[cnt,:,:,:]=deg_a
   
    
    deg_s=deg_c-deg_a;
    Dmat_sub[cnt,:,:,:]=deg_s
path_bf=data_path+ graphpath+'coh_degrees_15_resampled_2win_bothcon2_hubs_prob.mat'
ltc=scipy.io.loadmat(path_bf,mat_dtype=True); deg_bf=ltc['deg_final']
#for cc,bc in enumerate(trange):
#    b1=bc-1; b2=b1+20
#    Dmat_cncrt_r[:,:,:,cc]=np.subtract(Dmat_cncrt[:,:,:,b1:b2].mean(axis=3),Dmat_cncrt[:,:,:,20:40].mean(axis=3))
#    Dmat_abs_r[:,:,:,cc]=np.subtract(Dmat_abs[:,:,:,b1:b2].mean(axis=3),Dmat_abs[:,:,:,20:40].mean(axis=3))
#    Dmat_sub_r[:,:,:,cc]=np.subtract(Dmat_cncrt[:,:,:,b1:b2] .mean(axis=3),Dmat_abs[:,:,:,b1:b2] .mean(axis=3))
#print "resample fin"    
MTc=np.zeros(deg_bf.shape)
MTa=np.zeros(deg_bf.shape)
MTs=np.zeros(deg_bf.shape)
MPc=np.ones(deg_bf.shape)
MPa=np.ones(deg_bf.shape)
MPs=np.ones(deg_bf.shape)
for ii in range(MTc.shape[0]):
    print ii
    for jj in range(MTc.shape[1]):
        for kk in range(MTc.shape[2]):
            #if deg_bf[ii,jj,kk]>0:
#                b1=trange[kk]; b2=b1+20
#                dchunck=Dmat_cncrt[:,ii:ii+1,jj,b1:b2].mean(axis=2)-Dmat_abs[:,ii:ii+1,jj,b1:b2].mean(axis=2)#;dchunck=dchunck.reshape(dchunck.shape[0]*dchunck.shape[2],dchunck.shape[1])
                dchunck=Dmat_cncrt[:,ii:ii+1,jj,kk]-Dmat_abs[:,ii:ii+1,jj,kk]#;dchunck=dchunck.reshape(dchunck.shape[0]*dchunck.shape[2],dchunck.shape[1])

#                MTc[ii,jj,kk],MPc[ii,jj,kk],h=mne.stats.permutation_t_test(Dmat_cncrt_r[:,ii:ii+1,jj,kk], tail=0,n_permutations=1024)
#                MTa[ii,jj,kk],MPa[ii,jj,kk],h=mne.stats.permutation_t_test(Dmat_abs_r[:,ii:ii+1,jj,kk], tail=0,n_permutations=1024)
#                MTs[ii,jj,kk],MPs[ii,jj,kk],h=mne.stats.permutation_t_test(Dmat_sub_r[:,ii:ii+1,jj,kk], tail=0,n_permutations=1024)
                MTs[ii,jj,kk],MPs[ii,jj,kk],h=mne.stats.permutation_t_test(dchunck, tail=0,n_permutations=10000)
#dchunck=Dmat_cncrt[:,deg_bf.nonzero()[0],deg_bf.nonzero()[1],deg_bf.nonzero()[2]]-Dmat_abs[:,deg_bf.nonzero()[0],deg_bf.nonzero()[1],deg_bf.nonzero()[2]]#
#MTs,MPs,h=mne.stats.permutation_t_test(dchunck, tail=0,n_permutations=10000)
deg_bfn=deg_bf.reshape(deg_bf.shape[0],deg_bf.shape[1]*deg_bf.shape[2])
deg_bfns=np.sum(deg_bfn,axis=1)
MTs1=MTs.copy()
MTs[MPs>0.05]=0
MTsn=MTs.reshape(MTs.shape[0],MTs.shape[1]*MTs.shape[2])
MTsnp=MTsn.copy(); MTsnp[MTsnp<0]=0
MTsnn=MTsn.copy(); MTsnn[MTsnn>0]=0
MTsnsp=np.sum(MTsnp,axis=1)
MTsnsn=np.sum(MTsnn,axis=1)

#MCsn=range(deg_bfns.shape[0])
#MPsn=range(deg_bfns.shape[0])
#from mne.stats import permutation_t_test as stat_fun
#for ii in range(MTc.shape[0]):
#    print ii
#    if deg_bfns[ii]>0:
#        dchunck=Dmat_cncrt[:,ii,:,1::2]-Dmat_abs[:,ii,:,1::2]
#        Tobs,MCsn[ii],MPsn[ii],h=mne.stats.permutation_cluster_1samp_test(dchunck, tail=0,n_permutations=1024,out_type='indices',stat_fun=stat_fun)
##########
#Nper=1024        
##MTper=np.zeros((MTc.shape[0],len(range(Nper))))
#MTmaxP=np.zeros((len(range(Nper)),1))
#MTmaxN=np.zeros((len(range(Nper)),1))
#
#for per in range(Nper):
#    print per
#    MTsp=np.zeros(deg_bf.shape)
#    MPsp=np.ones(deg_bf.shape)
#    ranvec=0.25*np.random.randn(Dmat_cncrt.shape[0],1)
#    ranvec[ranvec>0]=1; ranvec[ranvec<0]=-1
#    for ii in range(MTc.shape[0]):
#        for jj in range(MTc.shape[1]):
#            for kk in range(MTc.shape[2]):
#                if deg_bf[ii,jj,kk]>0:
#                    
#    #                b1=trange[kk]; b2=b1+20
#    #                dchunck=Dmat_cncrt[:,ii:ii+1,jj,b1:b2].mean(axis=2)-Dmat_abs[:,ii:ii+1,jj,b1:b2].mean(axis=2)#;dchunck=dchunck.reshape(dchunck.shape[0]*dchunck.shape[2],dchunck.shape[1])
#                    dchunck=Dmat_cncrt[:,ii:ii+1,jj,kk]-Dmat_abs[:,ii:ii+1,jj,kk]#;dchunck=dchunck.reshape(dchunck.shape[0]*dchunck.shape[2],dchunck.shape[1])
#                    dchunck=dchunck*ranvec
#                    
#    #                MTc[ii,jj,kk],MPc[ii,jj,kk],h=mne.stats.permutation_t_test(Dmat_cncrt_r[:,ii:ii+1,jj,kk], tail=0,n_permutations=1024)
#    #                MTa[ii,jj,kk],MPa[ii,jj,kk],h=mne.stats.permutation_t_test(Dmat_abs_r[:,ii:ii+1,jj,kk], tail=0,n_permutations=1024)
#    #                MTs[ii,jj,kk],MPs[ii,jj,kk],h=mne.stats.permutation_t_test(Dmat_sub_r[:,ii:ii+1,jj,kk], tail=0,n_permutations=1024)
#                    MTsp[ii,jj,kk],MPsp[ii,jj,kk],h=mne.stats.permutation_t_test(dchunck, tail=0,n_permutations=10000)
#    MTsp[MPsp>0.05]=0
#    MTper=MTsp.reshape(MTsp.shape[0],MTsp.shape[1]*MTsp.shape[2])    
#    MTperp=MTper.copy(); MTperp[MTperp<0]=0
#    MTpern=MTper.copy(); MTpern[MTpern>0]=0
#    MTmaxP[per]=MTperp.max()#np.sum(MTperp,axis=1).max()
#    MTmaxN[per]=MTpern.min()#np.sum(MTpern,axis=1).min()
#path_MTn=data_path+ graphpath+'perN1024single_coh_degree_15_resampled_10win_50ms_bothcon2_hubs_prob.mat'
#scipy.io.savemat(path_MTn,{'MTmaxN':MTmaxN}) 
##ltc=scipy.io.loadmat(path_MTn,mat_dtype=True); MTmaxN=ltc['MTmaxN'] 
#path_MTp=data_path+ graphpath+'perP1024single_coh_degree_15_resampled_10win_50ms_bothcon2_hubs_prob.mat'
#scipy.io.savemat(path_MTp,{'MTmaxP':MTmaxP}) 
##ltc=scipy.io.loadmat(path_MTp,mat_dtype=True); MTmaxP=ltc['MTmaxP'] 
#sigsp=np.zeros(MTsnsp.shape)
#sigsn=np.zeros(MTsnsn.shape)
#for ii in range(len(sigsn)):
#    sigsp[ii]=np.sum(MTsnsp[ii]<MTmaxP)/len(MTmaxP)
#    sigsn[ii]=np.sum(MTsnsn[ii]>MTmaxN)/len(MTmaxP)
#
#sigsp=np.zeros(MTsnp.shape)
#sigsn=np.zeros(MTsnn.shape)
#for ii in range(sigsn.shape[0]):
#    for jj in range(sigsn.shape[1]):
#        for kk in range(sigsn.shape[2]):
#            sigsp[ii,jj,kk]=np.sum(MTsnp[ii,jj,kk]<MTmaxP)/len(MTmaxP)
#            sigsn[ii,jj,kk]=np.sum(MTsnn[ii,jj,kk]>MTmaxN)/len(MTmaxP)

                
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

#pvalc_c=MPc.copy()
#Pvalc=-10*np.log10(pvalc_c)
#Pvalc[Pvalc<-10*np.log10(0.05)]=0
#Pvalc[MTc<0]=-Pvalc[MTc<0]
#
#pvala_c=MPa.copy()
#Pvala=-10*np.log10(pvala_c)
#Pvala[Pvala<-10*np.log10(0.05)]=0
#Pvala[MTa<0]=-Pvala[MTa<0]

pvals_c=MPs.copy()
Pvals=-10*np.log10(pvals_c)
Pvals[Pvals<-10*np.log10(0.05)]=0
Pvals[MTs<0]=-Pvals[MTs<0]


      
fsave_vertices = [np.arange(10242), np.arange(10242)]
#Xc=np.zeros((20484, Pvalc.shape[2]))
#Xa=np.zeros((20484, Pvala.shape[2]))
Xs=np.zeros((20484, Pvals.shape[2]))


print "deg files loaded"
labelss = mne.read_labels_from_annot(subject='fsaverage', parc='rezvan_parc_mne_ctf_mn', subjects_dir=data_path)
label_names=[label.name for label in labelss]
labelss.remove(labelss[-2])
labelss.remove(labelss[-1])
for bn,band in zip(range(Pvals.shape[1]),['Theta','Alpha','Beta','Gamma']):
    print band
    for ccc,label in enumerate(labelss):  
        print ccc
        bb=stc_cncrt.in_label(label)
        if ccc%2==0:
            for cnl in bb.lh_vertno:
#                Xc[cnl,:]=Pvalc[ccc,bn,:]
#                Xa[cnl,:]=Pvala[ccc,bn,:]
                Xs[cnl,:]=Pvals[ccc,bn,:]
    
        else: 
            for cnr in bb.rh_vertno:
#                Xc[cnr+10242,:]=Pvalc[ccc,bn,:]
#                Xa[cnr+10242,:]=Pvala[ccc,bn,:]
                Xs[cnr+10242,:]=Pvals[ccc,bn,:]
                
        
       
    tmin1=5
    tstep1=5
    vertices_to = [np.arange(10242), np.arange(10242)]
#    Xc_stc = mne.SourceEstimate(Xc, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
#    Xa_stc = mne.SourceEstimate(Xa, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
    Xs_stc = mne.SourceEstimate(Xs, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
    
#    out_file_Xc=out_path + 'Perttest_10logPval_firstMorphed_SemLoc_icaclean_degree_proportional_25_pca_coh_concretebase_meanresampled_5_550_10ms_5overlap_'+band
#    Xc_stc.save(out_file_Xc)
#    
#    out_file_Xa=out_path + 'Perttest_10logPval_firstMorphed_SemLoc_icaclean_degree_proportional_25_pca_coh_abstractbase_meanresampled_5_550_10ms_5overlap_'+band
#    Xa_stc.save(out_file_Xa)
#    
    out_file_Xs=out_path  + 'Perttest_10logPval_firstMorphed_SemLoc_icaclean_hubs_proportional_15_pca_coh_subtract_meanresampled_50_550_100ms_0overlap_'+band
    Xs_stc.save(out_file_Xs)
# 
sigs=MPs[MPs<0.05]
sorted(sigs)<0.004*(np.arange(len(sigs))+1)
#out_path='/imaging/rf02/TypLexMEG/icaanalysis_results/jpg/time_resolved/ROI_connectivity/graph/coh_degforpca_hubs_meanz.mat';
#ltc=scipy.io.loadmat(out_path,mat_dtype=True); deg_pca=ltc['deg_pca']
#deg_pca=deg_pca[deg_bf>0,:];
#out_path='/imaging/rf02/TypLexMEG/icaanalysis_results/jpg/time_resolved/ROI_connectivity/graph/coh_degforpcaf_hubs_meanz.mat';
#scipy.io.savemat(out_path,{'deg_pca':deg_pca})