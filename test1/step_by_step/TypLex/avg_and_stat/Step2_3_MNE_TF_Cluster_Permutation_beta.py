"""
=========================================================
TF representation of TypLex for frequency bands in source labels
=========================================================

"""
# Authors: Alexandre Gramfort <gramfort@nmr.mgh.harvard.edu>
#
# License: BSD (3-clause)


print(__doc__)

import sys

sys.path.append('/imaging/local/software/python_packages/scikit-learn/v0.16.1')
sys.path.append('/imaging/local/software/mne_python/latest_v0.9')

# for qsub
# (add ! here if needed) 
sys.path.append('/imaging/local/software/anaconda/latest/x86_64/bin/python')
sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
###
sys.path.insert(1,'/imaging/local/software/python_packages/scikit-learn/v0.13.1/lib.linux-x86_64-2.7/sklearn/externals')
import joblib


import numpy as np
import mne
import os
import scipy
from mne import ( spatial_tris_connectivity, grade_to_tris)
from mne.stats import (spatio_temporal_cluster_1samp_test, summarize_clusters_stc)
import scipy.io as scio
###############################################################################
data_path = '/imaging/rf02/TypLexMEG/' # root directory for your MEG data
os.chdir(data_path)
subjects_dir = '/imaging/rf02/TypLexMEG/'    # where your MRI subdirectories are
# where event-files are
out_path = '/imaging/rf02/TypLexMEG/icaanalysis_results/typlex/stc/permutation/power/final2/' #
if not os.path.exists(out_path):
    os.makedirs(out_path)

label_path = '/imaging/rf02/TypLexMEG/fsaverage/label'

inv_path = '/imaging/rf02/TypLexMEG/'
inv_fname = 'InverseOperator_EMEG-inv.fif'

# get indices for subjects to be processed from command line input
# 
print sys.argv
p_inds = []
for ss in sys.argv[1:]:
   p_inds.append( int( ss ) )

subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16] # removed
#p_inds=[0,1,2,3,4,5,6,7,8,9,10]
p_list=[0.1,0.05,0.045,0.04,0.03,0.025,0.01,0.008,0.005,0.002,0.001,0.0005,0.0001]
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
for ss in p_inds:
 ll.append(p_list[ss])

print "ll:"
print ll



# Labels/ROIs
#labellist = ['atlleft-lh']#, 'atlleftt-rh'
tmin, tmax = -0.5, 0.7
reject = dict(eeg=120e-6, grad=200e-12, mag=4e-12)
event_ids = {'cncrt_wrd': 1, 'abs_wrd': 2}
# frequency bands with variable number of cycles for wavelets
stim_delay = 0.034 # delay in s
frequencies = np.arange(8, 45, 1)
n_cycles = frequencies / float(7)
n_cycles[frequencies<=20] = 2
n_cycles[frequencies<=20] += np.arange(0,1,0.077)
    # n_cycles[np.logical_and((frequencies>10), (frequencies<=20))] = 3
labellist = ['atlleft-lh', 'atlright-rh','medialtempright-rh','medialtempleft-lh']



ii=-1
X=np.zeros((20484,5,17))
Xc=np.zeros((20484,5,17))
Xa=np.zeros((20484,5,17))
#X_beta=np.zeros((20484,16))
b='beta'
for p_threshold in ll:
    for meg in list_all:
    	ii=ii+1;
    	print ii
    
    ### word
    	
    	fname_cncrt = data_path + meg + 'mineMorphed_ico_TypLex_ica_word_Power2_normori_ratio_equalized_m500_700_200hz_' +b#'mineMorphed_ico_TypLex_icaclean_word_Power_ratio_equalized_m500_700_' + b
    	stc_cncrt = mne.read_source_estimate(fname_cncrt)
    	#abslt_cncrt= np.absolute(stc_cncrt.copy().crop(.3,.4).data)
    	fname_abs = data_path + meg + 'mineMorphed_ico_TypLex_ica_nonword_Power2_normori_ratio_equalized_m500_700_200hz_' +b#'mineMorphed_ico_TypLex_icaclean_nonword_Power_ratio_equalized_m500_700_' + b
    	stc_abs = mne.read_source_estimate(fname_abs)
    	"""
    	fname_sub = data_path + meg + 'Morphed_TypLex_icaclean_Subtract_Power_ratio_equalized_m500_700_theta'
    	stc_subtract = mne.read_source_estimate(fname_sub)
    	#abslt_abs= np.absolute(stc_abs.copy().crop(.3,.4).data)	
    	#stc_subtract=np.subtract(stc_cncrt,stc_abs)
    	"""
    	b2=0.05
    	for cc in range(5):
    		b1=b2+0.0001; b2=b1+0.1-0.0001
    		print b1; print b2
    		Xc1=np.abs(stc_cncrt.copy().crop(b1,b2).mean().data.squeeze())
    		Xa1=np.abs(stc_abs.copy().crop(b1,b2).mean().data.squeeze())
    		X[:,cc,ii]=np.subtract(Xc1,Xa1)#np.subtract(np.abs(Xc1-1),np.abs(Xa1-1))#np.log(np.divide(Xc,Xa))#
    		Xc[:,cc,ii]=Xc1.copy()#np.log(np.divide(Xc,Xa))#
    		Xa[:,cc,ii]=Xa1.copy()#np.log(np.divide(Xc,Xa))#
      
    	#X[:,1,:]=X[:,0,:].copy()
    	#stc_subtract1=stc_subtract.copy().crop(.05,.6)
    	#stc_subtract1.resample(50)
    	#stc_cncrt1=stc_cncrt.copy().crop(.15,.351)
    	#stc_cncrt1.resample(50)
    	#X[:,:,ii]=stc_subtract1.data#.squeeze() for mean
    	#X[:,:,ii]=np.absolute(stc_subtract1.data)
      ####for PCA!!!
#    print X.shape
#    path_mtx=data_path+'TF_orig_words_'+b+'.mat'
#    scio.savemat(path_mtx,{'Xc':Xc})
#    path_mtx=data_path+'TF_orig_nonwords_'+b+'.mat'
#    scio.savemat(path_mtx,{'Xa':Xa})     
     
     
    	"""
    	X[:,0,ii]=bb.squeeze()
    	stc_subtract1=stc_subtract.crop(.25,.45).copy()
    	stc_subtract1=stc_subtract1.mean()
    	bb=stc_subtract1.data
    	X[:,1,ii]=bb.squeeze()
    	"""
    	#X_beta[:,ii]=stc_subtract.data[:,1]	
    	
    #X = np.transpose(X, [2, 1, 0])	
    #X_beta=X_beta.transpose()
    #X_beta=X_beta.transpose()	
    #t_beta,p_beta,H0_beta = clu_beta = mne.stats.permutation_t_test(X_beta, n_permutations=10000, tail=0, n_jobs=1, verbose=None)
    #t_beta,p_beta,H0_beta = clu_beta = mne.stats.permutation_t_test(X_beta, n_permutations=10000, tail=0, n_jobs=1, verbose=None)
    print('Computing connectivity.')
    
    fname_label = label_path + '/' + 'toremove_wbspokes-lh.label'; labelL = mne.read_label(fname_label)
    fname_label = label_path + '/' + 'toremove_wbspokes-rh.label'; labelR = mne.read_label(fname_label)
    labelss=labelL+labelR
    bb=stc_cncrt.in_label(labelss)
    fsave_vertices = [np.arange(10242), np.arange(10242)]
    nnl=np.in1d(fsave_vertices[0],bb.lh_vertno)
    nnr=np.in1d(fsave_vertices[1],bb.rh_vertno)
    spatial_exclude=np.hstack((fsave_vertices[0][nnl], fsave_vertices[0][nnr]+10242))
    
    connectivity = spatial_tris_connectivity(grade_to_tris(5))
    fsave_vertices = [np.arange(10242), np.arange(10242)]
    n_permutations=5000
    
    #    Note that X needs to be a multi-dimensional array of shape
    #    samples (subjects) x time x space, so we permute dimensions
    X = np.transpose(X, [2, 1, 0])
    
    #    Now let's actually do the clustering. This can take a long time...
    #    Here we set the threshold quite high to reduce computation.
    #p_threshold = 0.05
    n_subjects=17
    t_threshold = -scipy.stats.distributions.t.ppf(p_threshold/2., n_subjects - 1)
    print('Clustering.')
    max_step=1
    T_obs, clusters, cluster_p_values, H0 = clu = spatio_temporal_cluster_1samp_test(X, connectivity=connectivity, n_jobs=4, threshold=t_threshold,n_permutations=n_permutations,tail=0,t_power=1, step_down_p=0.05, spatial_exclude=spatial_exclude)#, max_step=5)
    
    #    Now select the clusters that are sig. at p < 0.05 (note that this value
    #    is multiple-comparisons corrected).

#    good_cluster_inds = np.where(cluster_p_values < 0.05)[0]
#    print cluster_p_values.min()
    tstep2=0.1
    tmin2=0.05
    clus_p_values=(p_threshold/0.05)*cluster_p_values
               
    good_cluster_inds = np.where(cluster_p_values<=0.05)[0]#np.intersect1d(np.where(clus_p_values <= 0.01)[0], np.where(cluster_p_values<=0.1)[0])
    print cluster_p_values.min(), clus_p_values.min()
        
    fsave_vertices = [np.arange(10242), np.arange(10242)]
    if good_cluster_inds.size>0:
        p_thresh=np.max(cluster_p_values[good_cluster_inds])#cluster_p_values.min()+0.000001
        	
    if cluster_p_values.min()<=0.1:#clus_p_values.min()<=0.01 and cluster_p_values.min()<=0.1:
        print p_thresh
                
        stc_all_cluster_vis = summarize_clusters_stc(clu, tstep=tstep2, p_thresh=p_thresh+0.000001, vertno=fsave_vertices, subject='fsaverage')
        #t_stc = mne.SourceEstimate(T_obs, vertices=fsave_vertices, tmin=tmin2, tstep=tstep2, subject='average', verbose=None)
        
        out_file1=out_path + 'Per_clusp'+str(p_threshold)[2:]+'_p'+str(cluster_p_values.min())[2:]+'_TL_ica_Subtract_Power_ratio_ico_normori_eq_50_550_100ms_sx_ms' + str(max_step) +b
        stc_all_cluster_vis.save(out_file1)
        #out_file2=data_path + 'Permutation_TypLex_Subtract_Coherence_Beta_150_350_' + label_name[0:-3]
        #clu_beta.save(out_file2)
        
        
        Matx=np.zeros((20484,5))
        T=np.divide(T_obs,np.absolute(T_obs))
        for cc in range(good_cluster_inds.shape[0]):
        	Matx[clusters[good_cluster_inds[cc]][1],clusters[good_cluster_inds[cc]][0]]=T[clusters[good_cluster_inds[cc]][0],clusters[good_cluster_inds[cc]][1]]
        
        tmin1=50
        tstep1=100
        vertices_to = [np.arange(10242), np.arange(10242)]
        matx_stc = mne.SourceEstimate(Matx, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
        out_file2=out_path + 'Per_sw_clusp'+str(p_threshold)[2:]+'_p'+str(p_thresh)[2:]+'_TL_ica_Subtract_Power_ratio_normori_eq_50_550_100ms_sx_ms' + str(max_step) +b
        matx_stc.save(out_file2)

