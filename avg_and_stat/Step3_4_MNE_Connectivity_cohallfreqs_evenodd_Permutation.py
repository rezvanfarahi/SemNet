# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 16:07:07 2016

@author: rf02
"""


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
#sys.path.append('/imaging/local/software/python_packages/nibabel/1.3.0')
sys.path.insert(1,'/imaging/rf02/scikit-learn-0.16.1')
#sys.path.append('/imaging/local/software/python_packages/scikit-learn/v0.14.1')
#sys.path.insert(1,'/imaging/rf02/Semnet/mne_python_0.14')
sys.path.insert(1,'/imaging/local/software/mne_python/latest_v0.11')

# for qsub
# (add ! here if needed) 
#sys.path.insert(1,'/imaging/local/software/anaconda/1.9.1/x86_64/envs/py3k')#bin/python')
sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
sys.path.insert(1,'/home/rf02/.local/lib/python2.7/site-packages')
#from sklearn.decomposition import PCA
####
sys.path.insert(1,'/imaging/local/software/python_packages/scikit-learn/v0.13.1/lib.linux-x86_64-2.7/sklearn/externals')
import joblib
import numpy as np
import mne
#reload(mne)
import os
import scipy.io
import scipy.stats as scist

from mne import ( spatial_tris_connectivity, grade_to_tris)
from mne.stats import (spatio_temporal_cluster_1samp_test,
                       summarize_clusters_stc)


###############################################################################
data_path = '/imaging/rf02/Semnet/' # root directory for your MEG data
out_path = '/imaging/rf02/Semnet/stc/permutation/connectivity/VisualHand/maskedlh/top10verts_separatetimes/'    # where event files are
if not os.path.exists(out_path):
    os.makedirs(out_path)
label_path1 = '/imaging/rf02/TypLexMEG/fsaverage/label'

# 
print sys.argv
p_inds = []
for ss in sys.argv[1:]:
   p_inds.append( int( ss ) )

subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18] # removed
#p_inds=[0,1,2,3,4,5,6,7,8,9,10]
#p_list=[0.05,0.025,0.01,0.008,0.005,0.002,0.001]#[0.005,0.01,0.025, 0.05]#,0.045,0.04,0.03,0.025,0.02,0.01,0.008,0.005,0.002,0.001,0.0005,0.0001]
p_list=[0.1,0.05,0.045,0.04,0.03,0.025,0.02,0.01,0.008,0.005,0.002,0.001, 0.0008,0.0005,0.0001]#[0.05,0.025,0.01,0.008,0.005,0.002,0.001]#[0.005,0.01,0.025, 0.05]#,0.045,0.04,0.03,0.025,0.02,0.01,0.008,0.005,0.002,0.001,0.0005,0.0001]

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
            '/meg16_0122/160707/', #22 
            '/meg16_0125/160712/', #24 
            ]



# subjects names used for MRI data
subjects=  ['MRI_meg16_0030' ,#0
            'MRI_meg16_0032' ,#1
            'MRI_meg16_0034' ,#2
            'MRI_meg16_0035' ,#3
            'MRI_meg16_0042' ,#4
            'MRI_meg16_0045' ,#5
            'MRI_meg16_0052' ,#6
            'MRI_meg16_0056' ,#7
    	       'MRI_meg16_0069' ,#8
 	       'MRI_meg16_0070' ,#9
            'MRI_meg16_0072' ,#10
            'MRI_meg16_0073' ,#11
            'MRI_meg16_0075' ,#12
            'MRI_meg16_0078' ,#13
            'MRI_meg16_0082' ,#14
            'MRI_meg16_0086' ,#15
            'MRI_meg16_0097' ,#16
            'MRI_meg16_0122' ,#17
            'MRI_meg16_0125' ,#18
            ]

ll = []
for ss in p_inds:
 ll.append(p_list[ss])

print "ll:"
print ll



# Labels/ROIs
#labellist = ['atlleft-lh']#, 'atlleftt-rh'
tmin, tmax = -0.5, 0.7
reject = dict(eeg=120e-6, grad=200e-12, mag=4e-12) #eog=150e-6, 
event_ids = {'cncrt_wrd': 1, 'abs_wrd': 2}
# frequency bands with variable number of cycles for wavelets
stim_delay = 0.034 # delay in s
frequencies = np.arange(8, 45, 1)
n_cycles = frequencies / float(7)
n_cycles[frequencies<=20] = 2
n_cycles[frequencies<=20] += np.arange(0,1,0.077)
    # n_cycles[np.logical_and((frequencies>10), (frequencies<=20))] = 3
#'lh.lateralorbitofrontal','rh.lateralorbitofrontal',
#labellist=['lhATL_latmed-lh','rhATL_latmed-rh','lhAG_hsmetafinal5-lh','lhposterior_inferiorfrontal-lh','lhmiddle_temporal-lh']#labellist = ['lh.ATL_latmed','lh.supramarginal','lh.posterior_inferiorfrontal','lh.middle_temporal']#labellist = [ 'lh.posterior_fusiform','rh.posterior_fusiform', 'lh.lateraloccipital','rh.lateraloccipital','lh.precentral','rh.precentral']#, 'lh.ATLsmall','rh.ATLsmall']#labellist = [ 'leftATL-lh','rightATL-rh']#, 'medialtempright-rh','medialtempleft-lh']'lh.postcentral','rh.postcentral',
#labellist2=['lhATL-lh','rhATL-rh','lhAG-lh','lhpIFG-lh','lhMTL-lh']#labellist = ['lh.ATL_latmed','lh.supramarginal','lh.posterior_inferiorfrontal','lh.middle_temporal']#labellist = [ 'lh.posterior_fusiform','rh.posterior_fusiform', 'lh.lateraloccipital','rh.lateraloccipital','lh.precentral','rh.precentral']#, 'lh.ATLsmall','rh.ATLsmall']#labellist = [ 'leftATL-lh','rightATL-rh']#, 'medialtempright-rh','medialtempleft-lh']'lh.postcentral','rh.postcentral',
#labellist1=['lh0_ATL_latmed_ctf', 'lh1_AG_hsmetafinal5_ctf', 'lh2_col_shape_c5_ctf', 'rh3_col_shape_c5_ctf', 'lh4_auditory_c_ctf', 'rh5_auditory_c_ctf', 'lh6_hand_c_ctf', 'rh7_hand_c_ctf']#'lh6_hand_c', 'rh7_hand_c']#'lh4_auditory_c', 'rh5_auditory_c']#, 'hand_c-lh', 'hand_c-rh']
labellist1 = ['lh0_ATL_ventromedial_top10verts_ctf','lh1_ATL_dorsolateral_top10verts_ctf','lh2_AG_anterior_top10verts_ctf','lh3_AG_posterior_top10verts_ctf','lh4_hand_postcentral_top10verts_ctf','lh5_hand_precentral_top10verts_ctf','lh6_col_shape_c5_top10verts_ctf']#,label_path+'5_auditory_c_ctf-rh.label',label_path+'6_hand_c_ctf-lh.label',label_path+'7_hand_c_ctf-rh.label']#,label_path+'2_col_shape_c5_ctf-lh.label',label_path+'3_col_shape_c5_ctf-rh.label',label_path+'4_auditory_c_ctf-lh.label',label_path+'5_auditory_c_ctf-rh.label',label_path+'6_hand_c_ctf-lh.label',label_path+'7_hand_c_ctf-rh.label']#[ 'lh.posterior_fusiform','rh.posterior_fusiform', 'lh.lateraloccipital','rh.lateraloccipital','lh.precentral','rh.precentral', 'lh.ATL_latmed','rh.ATL_latmed']#labellist = [ 'leftATL-lh','rightATL-rh']#, 'medialtempright-rh','medialtempleft-lh']'lh.postcentral','rh.postcentral',


#labellist_path = [label_path+'0_ATL_latmed_ctf-lh.label',label_path+'1_AG_hsmetafinal5_ctf-lh.label',label_path+'2_col_shape_c5_ctf-lh.label',label_path+'3_col_shape_c5_ctf-rh.label',label_path+'4_auditory_c_ctf-lh.label',label_path+'5_auditory_c_ctf-rh.label',label_path+'6_hand_c_ctf-lh.label',label_path+'7_hand_c_ctf-rh.label']#[ 'lh.posterior_fusiform','rh.posterior_fusiform', 'lh.lateraloccipital','rh.lateraloccipital','lh.precentral','rh.precentral', 'lh.ATL_latmed','rh.ATL_latmed']#labellist = [ 'leftATL-lh','rightATL-rh']#, 'medialtempright-rh','medialtempleft-lh']'lh.postcentral','rh.postcentral',

#labellist21=['lhATL-lh','lhAG-lh','lhvis-lh','rhvis-rh','lhaud-lh','rhaud-rh','lhhnd-lh','rhhnd-rh']
labellist21=['lhATL_ventromedial-lh','lhATL_dorsolateral-lh','lhAG_anterior-lh','lhAG_posterior-lh','lhhnd_postc-lh','lhhnd_prec-lh','lhvis-lh']#['lhATL','lhAG','lhIFG','lhMTG','lhvis','rhvis','lhaud','rhaud','lhhnd','rhhnd']

#labellist = [ 'lh.posterior_fusiform','rh.posterior_fusiform', 'lh.lateraloccipital','rh.lateraloccipital','lh.precentral','rh.precentral']#labellist = [ 'posterior_fusiform-lh','rh.lateralorbitofron-rh', 'rh.medialorbitofron-rh']#labellist = ['leftATL-rh','rightATL-lh']#, 'rightATL-rh']#'atlleft-lh']#, 'atlright-rh']#,'medialtempright-rh','medialtempleft-lh']
labellist=[labellist1[ii] for ii in np.array([0,1,2,3,4,5,6])]
labellist2=[labellist21[ii] for ii in np.array([0,1,2,3,4,5,6])]
event1='Visual'
event2='Hand'
n_subjects=len(list_all)
ntimes=4
nbands=4
bands=['Theta','Alpha','Beta','Gamma']
for p_threshold in ll:
    for label_name,label_name_sv in zip(labellist,labellist2):
            
            	ii=-1;
            	#X=np.zeros((20484,2,17))
#            	X=np.zeros((20484,2,17))
            	Xce=np.zeros((20484,ntimes,nbands,n_subjects))
            	Xae=np.zeros((20484,ntimes,nbands,n_subjects))
             
            	
            	for meg in list_all:
                     ii=ii+1;
                     print ii
                     for bcnt,band in enumerate(bands):
            		
            	
            ### Concrete
                		fname_event1e = data_path + meg + 'Morphed_ico_SemDec_evenodd_ica_equalized_'+event1+'_Mlttap_top10verts_ctf_Coh_'+band+'_150_450_200ms_' + label_name 
                		stc_event1e = mne.read_source_estimate(fname_event1e)#; stc_cncrt.data[stc_cncrt.data<0]=0 
                		
                		fname_event2e = data_path + meg + 'Morphed_ico_SemDec_evenodd_ica_equalized_'+event2+'_Mlttap_top10verts_ctf_Coh_'+band+'_150_450_200ms_' + label_name   
                		stc_event2e = mne.read_source_estimate(fname_event2e)#; stc_abs.data[stc_abs.data<0]=0 	
                		
                		data_event1e = np.abs(stc_event1e.data)#; data_cncrt_wpli2d = np.abs(stc_cncrt_wpli2d.data); data_cncrt_coh = np.abs(stc_cncrt_coh.data)
                		data_event2e = np.abs(stc_event2e.data)#; data_abs_wpli2d = np.abs(stc_abs_wpli2d.data); data_abs_coh = np.abs(stc_abs_coh.data)
                		Xce[:,:,bcnt,ii]=data_event1e[:,:ntimes]#np.dstack((data_cncrt_mi,data_cncrt_wpli2d,data_cncrt_coh))
                		Xae[:,:,bcnt,ii]=data_event2e[:,:ntimes]#np.dstack((data_abs_mi,data_abs_wpli2d,data_abs_coh))
                      
#            	Xcfe=np.hstack((np.squeeze(Xce[:,:,0,:]),np.squeeze(Xce[:,:,1,:]),np.squeeze(Xce[:,:,2,:]),np.squeeze(Xce[:,:,3,:]))) #mi,  coh
#            	Xafe=np.hstack((np.squeeze(Xae[:,:,0,:]),np.squeeze(Xae[:,:,1,:]),np.squeeze(Xae[:,:,2,:]),np.squeeze(Xae[:,:,3,:])))
            	Xcfe=np.hstack((np.squeeze(Xce[:,0,:,:]),np.squeeze(Xce[:,1,:,:]),np.squeeze(Xce[:,2,:,:]),np.squeeze(Xce[:,3,:,:]))) #mi,  coh
            	Xafe=np.hstack((np.squeeze(Xae[:,0,:,:]),np.squeeze(Xae[:,1,:,:]),np.squeeze(Xae[:,2,:,:]),np.squeeze(Xae[:,3,:,:])))
                   
            	X = np.subtract(Xcfe,Xafe)#np.subtract(Xcpca,Xapca)
            	X = np.transpose(X, [2, 1, 0])	
            	X = X[:,np.array([8,9,10,11]),:]#,7,6,5,4,8,9,10,11	
             
            	connectivity = spatial_tris_connectivity(grade_to_tris(5))
            	#p_threshold = 0.045
            	print p_threshold
            	n_subjects=n_subjects
            	t_threshold = -scist.distributions.t.ppf(p_threshold/2., n_subjects - 1)
            	fsave_vertices = [np.arange(10242), np.arange(10242)]
            	n_permutations=5000
            	# I'm replicating the columns because it makes the script faster. Have checked and it doesn't introduce 
            	
            	print('Computing connectivity.')
            	fname_label = label_path1 + '/' + 'toremove_wbspokes-lh.label'; labelL = mne.read_label(fname_label)                 
#            	fname_label = label_path + '/' + 'toremove_wbspokes-rh.label'; labelR = mne.read_label(fname_label)
            	label_path='/imaging/rf02/Semnet/masks'
#                                     
#            	fname_label = label_path + '/' + 'mask_vishand_lang_L-lh.label'; labelL = mne.read_label(fname_label)
            	fname_label = label_path + '/' + 'mask_vishand_lang_R-rh.label'; labelR = mne.read_label(fname_label)
            	labelss=labelL+labelR
            	bb=stc_event1e.in_label(labelss)
            	fsave_vertices = [np.arange(10242), np.arange(10242)]
            	nnl=np.in1d(fsave_vertices[0],bb.lh_vertno)
            	nnr=np.in1d(fsave_vertices[1],bb.rh_vertno)
            	spatial_exclude=np.hstack((fsave_vertices[0][nnl], fsave_vertices[0][nnr]+10242))
                
            
            	
            	print('Clustering theta'); max_step=1;
            	T_obst, clusterst, cluster_p_valuest, H0t = clut = spatio_temporal_cluster_1samp_test(X, connectivity=connectivity, n_jobs=4, n_permutations=n_permutations,  threshold=t_threshold, tail=0,t_power=1, spatial_exclude=spatial_exclude,step_down_p=0.05)#, max_step=1)
            	#    Now select the clusters that are sig. at p < 0.05 (note that this value
            	#    is multiple-comparisons corrected).
            	clus_p_values=(p_threshold/0.05)*cluster_p_valuest
               
            	good_cluster_indst =  np.where(cluster_p_valuest <= 0.1)[0]#np.where(clus_p_values <= 0.003125)[0]
            	print cluster_p_valuest.min(), clus_p_values.min()
            	
            	
            	
            
            	tstep=1e-3
            	fsave_vertices = [np.arange(10242), np.arange(10242)]
            	if good_cluster_indst.size>0:
                	p_thresh=np.max(cluster_p_valuest[good_cluster_indst])#cluster_p_valuest.min()+0.000001
            	
            	if cluster_p_valuest.min()<=0.1:
                    p_thresh=cluster_p_valuest.min()#0.1#0.05; print p_thresh
                    
                    stc_all_cluster_vist = summarize_clusters_stc(clut, tstep=tstep, p_thresh=p_thresh+0.000001, vertices=fsave_vertices, subject='fsaverage')
                    out_file1=out_path + 'ClusP_Morphed_ico_SD_ica_evenodd_Sub_ctf_ST_allbands_250_450_200ms_ms' + str(max_step) +'_'+ str(p_threshold)[2:5] + '_sx_' + str(p_thresh)[2:5]+ '_' + label_name_sv[0:-3]
                    stc_all_cluster_vist.save(out_file1)
                      	
                    Matx=np.zeros((20484,X.shape[1]))
                    T=np.divide(T_obst,np.absolute(T_obst))
                    for cc in range(good_cluster_indst.shape[0]):
                        Matx[clusterst[good_cluster_indst[cc]][1],clusterst[good_cluster_indst[cc]][0]]=T[clusterst[good_cluster_indst[cc]][0],clusterst[good_cluster_indst[cc]][1]]
                    
                    tmin1=0
                    tstep1=100
                    vertices_to = [np.arange(10242), np.arange(10242)]
                    matx_stc = mne.SourceEstimate(Matx, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
                    out_file2=out_path + 'ClusP_sw_Morphed_ico_SD_ica_evenodd_Sub_ctf_ST_allbands_250_450_200ms_ms' + str(max_step) +'_sx_'+ str(p_threshold)[2:5] + '_' + str(p_thresh)[2:5]+ '_' + label_name_sv[0:-3]
                    matx_stc.save(out_file2)
            	
