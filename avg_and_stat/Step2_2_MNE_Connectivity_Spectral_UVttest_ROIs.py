
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
sys.path.insert(1,'/imaging/local/software/mne_python/latest')
# for qsub
# (add ! here if needed) 
#sys.path.insert(1,'/imaging/local/software/anaconda/1.9.1/x86_64/envs/py3k')#bin/python')
sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
###
import numpy as np
import mne
from mne import io

import os
import scipy.io
from scipy import stats

from mne import (io, spatial_tris_connectivity, compute_morph_matrix,
                 grade_to_tris)
from mne.stats import (spatio_temporal_cluster_1samp_test,
                       summarize_clusters_stc)


###############################################################################
data_path =  '/imaging/rf02/Semnet/' # root directory for your MEG data
os.chdir(data_path)
subjects_dir =  '/imaging/rf02/Semnet/'    # where your MRI subdirectories are
# where event-files are
out_path = '/imaging/rf02/Semnet/stc/uvttest/connectivity/evenodd/'    # where event files are
if not os.path.exists(out_path):
    os.makedirs(out_path)

avg_path = '/imaging/rf02/Semnet/stc/GrandAverage/connectivity/evenodd/'    # where event files are
if not os.path.exists(avg_path):
    os.makedirs(avg_path)
label_path='/imaging/rf02/Semnet/stc/localisers/'
#inv_path = '/imaging/rf02/TypLexMEG/'
#inv_fname = 'InverseOperator_EMEG-inv.fif'

# get indices for subjects to be processed from command line input
# 

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


#original ROIs
#labellist_path = [label_path+'col_shape_c5-lh.label',label_path+'col_shape_c5-rh.label',label_path+'auditory_c-lh.label',label_path+'auditory_c-rh.label',label_path+'hand_c-lh.label',label_path+'hand_c-rh.label']#[ 'lh.posterior_fusiform','rh.posterior_fusiform', 'lh.lateraloccipital','rh.lateraloccipital','lh.precentral','rh.precentral', 'lh.ATL_latmed','rh.ATL_latmed']#labellist = [ 'leftATL-lh','rightATL-rh']#, 'medialtempright-rh','medialtempleft-lh']'lh.postcentral','rh.postcentral',
##labellist_path = [label_path+'col_shape_c5-lh.label',label_path+'col_shape_c5-rh.label',label_path+'auditory_c-lh.label',label_path+'auditory_c-rh.label',label_path+'hand_c-lh.label',label_path+'hand_c-rh.label']#[ 'lh.posterior_fusiform','rh.posterior_fusiform', 'lh.lateraloccipital','rh.lateraloccipital','lh.precentral','rh.precentral', 'lh.ATL_latmed','rh.ATL_latmed']#labellist = [ 'leftATL-lh','rightATL-rh']#, 'medialtempright-rh','medialtempleft-lh']'lh.postcentral','rh.postcentral',
#labellist_spokes=[mne.read_label(label) for label in labellist_path]
#for thisl in labellist_spokes:
#    thisl.values.fill(1.0)
#
#fsaverage_path='/imaging/rf02/TypLexMEG'
#labelss = mne.read_labels_from_annot(subject='fsaverage', parc='rez_aparc_19oct16', subjects_dir=fsaverage_path)
#
#for lbl in labelss:
#    lbl.values.fill(1.0)
#label_names=[label.name for label in labelss]
#label_index=np.array([label_names.index('ATL_latmed-lh'),label_names.index('AG_hsmetafinal5-lh'),label_names.index('posterior_inferiorfrontal-lh'),label_names.index('middle_temporal-lh')])
#
#labellist_hubs=[labelss[ii] for ii in label_index]
#halfmax verts
#labellist_path = [label_path+'0_ATL_latmed_halfmaxverts_ctf-lh.label',label_path+'1_AG_hsmetafinal5_halfmaxverts_ctf-lh.label',label_path+'2_posterior_inferiorfrontal_halfmaxverts_ctf-lh.label',label_path+'3_middle_temporal_halfmaxverts_ctf-lh.label',label_path+'4_col_shape_c5_halfmaxverts_ctf-lh.label',label_path+'5_col_shape_c5_halfmaxverts_ctf-rh.label',label_path+'6_auditory_c_halfmaxverts_ctf-lh.label',label_path+'7_auditory_c_halfmaxverts_ctf-rh.label',label_path+'8_hand_c_halfmaxverts_ctf-lh.label',label_path+'9_hand_c_halfmaxverts_ctf-rh.label']#[ 'lh.posterior_fusiform','rh.posterior_fusiform', 'lh.lateraloccipital','rh.lateraloccipital','lh.precentral','rh.precentral', 'lh.ATL_latmed','rh.ATL_latmed']#labellist = [ 'leftATL-lh','rightATL-rh']#, 'medialtempright-rh','medialtempleft-lh']'lh.postcentral','rh.postcentral',
#labellist_path = [label_path+'0_ATL_ventromedial_top10verts_ctf-lh.label',label_path+'1_ATL_dorsolateral_top10verts_ctf-lh.label',label_path+'2_AG_anterior_top10verts_ctf-lh.label',label_path+'3_AG_posterior_top10verts_ctf-lh.label',label_path+'4_hand_postcentral_top10verts_ctf-lh.label',label_path+'5_hand_precentral_top10verts_ctf-lh.label',label_path+'6_col_shape_c5_top10verts_ctf-lh.label']#,label_path+'5_auditory_c_ctf-rh.label',label_path+'6_hand_c_ctf-lh.label',label_path+'7_hand_c_ctf-rh.label']#,label_path+'2_col_shape_c5_ctf-lh.label',label_path+'3_col_shape_c5_ctf-rh.label',label_path+'4_auditory_c_ctf-lh.label',label_path+'5_auditory_c_ctf-rh.label',label_path+'6_hand_c_ctf-lh.label',label_path+'7_hand_c_ctf-rh.label']#[ 'lh.posterior_fusiform','rh.posterior_fusiform', 'lh.lateraloccipital','rh.lateraloccipital','lh.precentral','rh.precentral', 'lh.ATL_latmed','rh.ATL_latmed']#labellist = [ 'leftATL-lh','rightATL-rh']#, 'medialtempright-rh','medialtempleft-lh']'lh.postcentral','rh.postcentral',
labellist_path = [label_path+'0_ATL_ventromedial_halfmaxverts_ctf-lh.label',label_path+'1_ATL_dorsolateral_halfmaxverts_ctf-lh.label',label_path+'2_AG_anterior_halfmaxverts_ctf-lh.label',label_path+'3_AG_posterior_halfmaxverts_ctf-lh.label',label_path+'4_hand_postcentral_halfmaxverts_ctf-lh.label',label_path+'5_hand_precentral_halfmaxverts_ctf-lh.label',label_path+'8_hand_c_halfmaxverts_ctf-lh.label']#,label_path+'5_auditory_c_ctf-rh.label',label_path+'6_hand_c_ctf-lh.label',label_path+'7_hand_c_ctf-rh.label']#,label_path+'2_col_shape_c5_ctf-lh.label',label_path+'3_col_shape_c5_ctf-rh.label',label_path+'4_auditory_c_ctf-lh.label',label_path+'5_auditory_c_ctf-rh.label',label_path+'6_hand_c_ctf-lh.label',label_path+'7_hand_c_ctf-rh.label']#[ 'lh.posterior_fusiform','rh.posterior_fusiform', 'lh.lateraloccipital','rh.lateraloccipital','lh.precentral','rh.precentral', 'lh.ATL_latmed','rh.ATL_latmed']#labellist = [ 'leftATL-lh','rightATL-rh']#, 'medialtempright-rh','medialtempleft-lh']'lh.postcentral','rh.postcentral',

labellist_spokes=[mne.read_label(label) for label in labellist_path]
for thisl in labellist_spokes:
    thisl.values.fill(1.0)

labellist_target_pre=labellist_spokes
#labellist_sv_target_pre=['lhATL','lhAG','lhIFG','lhMTG','lhvis','rhvis','lhaud','rhaud','lhhnd','rhhnd']
#labellist_sv_target_pre=['lhATL','lhAG','lhIFG','lhMTG','lhvis','rhvis','lhaud','rhaud','lhhnd','rhhnd']
labellist_sv_target_pre=['lhATL_ventromedial','lhATL_dorsolateral','lhAG_anterior','lhAG_posterior','lhhnd_postc','lhhnd_prec','lhvis']
#'lh.lateralorbitofrontal','rh.lateralorbitofrontal']#,
#labellist=['lhATL_latmed-lh','rhATL_latmed-rh','lhAG_hsmetafinal5-lh','lhposterior_inferiorfrontal-lh','lhmiddle_temporal-lh']#labellist = ['lh.ATL_latmed','lh.supramarginal','lh.posterior_inferiorfrontal','lh.middle_temporal']#labellist = [ 'lh.posterior_fusiform','rh.posterior_fusiform', 'lh.lateraloccipital','rh.lateraloccipital','lh.precentral','rh.precentral', 'lh.ATLsmall','rh.ATLsmall']#labellist = [ 'leftATL-lh','rightATL-rh']#, 'medialtempright-rh','medialtempleft-lh']'lh.postcentral','rh.postcentral',
#label_path+'lh0_posterior_inferiorfrontal_ctf-lh.label',label_path+'lh1_middle_temporal_ctf-lh.label'
#labellist_path = [label_path+'0_ATL_latmed_ctf-lh.label',label_path+'1_AG_hsmetafinal5_ctf-lh.label',label_path+'2_col_shape_c5_ctf-lh.label',label_path+'3_col_shape_c5_ctf-rh.label',label_path+'4_auditory_c_ctf-lh.label',label_path+'5_auditory_c_ctf-rh.label',label_path+'6_hand_c_ctf-lh.label',label_path+'7_hand_c_ctf-rh.label']#[ 'lh.posterior_fusiform','rh.posterior_fusiform', 'lh.lateraloccipital','rh.lateraloccipital','lh.precentral','rh.precentral', 'lh.ATL_latmed','rh.ATL_latmed']#labellist = [ 'leftATL-lh','rightATL-rh']#, 'medialtempright-rh','medialtempleft-lh']'lh.postcentral','rh.postcentral',

#labellist1=['lh0_ATL_latmed_ctf', 'lh1_AG_hsmetafinal5_ctf','lh0_posterior_inferiorfrontal_ctf', 'lh1_middle_temporal_ctf', 'lh2_col_shape_c5_ctf', 'rh3_col_shape_c5_ctf', 'lh4_auditory_c_ctf', 'rh5_auditory_c_ctf', 'lh6_hand_c_ctf', 'rh7_hand_c_ctf']#'lh6_hand_c', 'rh7_hand_c']#'lh4_auditory_c', 'rh5_auditory_c']#, 'hand_c-lh', 'hand_c-rh']
labellist1 = ['lh0_ATL_ventromedial_top10verts_ctf','lh1_ATL_dorsolateral_top10verts_ctf','lh2_AG_anterior_top10verts_ctf','lh3_AG_posterior_top10verts_ctf','lh4_hand_postcentral_top10verts_ctf','lh5_hand_precentral_top10verts_ctf','lh6_col_shape_c5_top10verts_ctf']#,label_path+'5_auditory_c_ctf-rh.label',label_path+'6_hand_c_ctf-lh.label',label_path+'7_hand_c_ctf-rh.label']#,label_path+'2_col_shape_c5_ctf-lh.label',label_path+'3_col_shape_c5_ctf-rh.label',label_path+'4_auditory_c_ctf-lh.label',label_path+'5_auditory_c_ctf-rh.label',label_path+'6_hand_c_ctf-lh.label',label_path+'7_hand_c_ctf-rh.label']#[ 'lh.posterior_fusiform','rh.posterior_fusiform', 'lh.lateraloccipital','rh.lateraloccipital','lh.precentral','rh.precentral', 'lh.ATL_latmed','rh.ATL_latmed']#labellist = [ 'leftATL-lh','rightATL-rh']#, 'medialtempright-rh','medialtempleft-lh']'lh.postcentral','rh.postcentral',

labellist_sv1=['lhATL_ventromedial','lhATL_dorsolateral','lhAG_anterior','lhAG_posterior','lhhnd_postc','lhhnd_prec','lhvis']#['lhATL','lhAG','lhIFG','lhMTG','lhvis','rhvis','lhaud','rhaud','lhhnd','rhhnd']

#labellist = [ 'lh.posterior_fusiform','rh.posterior_fusiform']#,'lh.precentral','rh.precentral', 'lh.lateraloccipital','rh.lateraloccipital']#labellist = [ 'posterior_fusiform-lh','rh.lateralorbitofron-rh', 'rh.medialorbitofron-rh']#labellist = ['leftATL-rh','rightATL-lh']#, 'rightATL-rh']#'atlleft-lh']#, 'atlright-rh']#,'medialtempright-rh','medialtempleft-lh']
conn_methods=['Coherence']#,'DWPLI']
bands=['Theta','Alpha','Beta','Gamma']
ntimes=4
pick_seeds=np.array([6])
pick_targets=np.array([0,1,2,3])
labellist_seed=[labellist1[ii] for ii in pick_seeds]
labellist_target=[labellist_target_pre[ii] for ii in pick_targets]
labellist_sv_seed=[labellist_sv1[ii] for ii in pick_seeds]
labellist_sv_target=[labellist_sv_target_pre[ii] for ii in pick_targets]
event1='Visual'
event2='Hear'
for seed_name, seed_name_sv in zip(labellist_seed, labellist_sv_seed):
    print seed_name
    X=np.zeros((19,ntimes,len(labellist_target),len(bands)))  
    MT=np.zeros((1,ntimes,len(labellist_target),len(bands)))   
    tcnt=-1                     
    for target_label, target_name in zip(labellist_target, labellist_sv_target):
        tcnt=tcnt+1
        for bandc,band in enumerate(bands):
            
            ii=-1; print band
            
#            Xa=np.zeros((20484,3,19))
            #    	X_theta=np.zeros((20484,17))
            #    	X_alpha=np.zeros((20484,17))
            #    	X_beta=np.zeros((20484,17))
            #    	X_gamma=np.zeros((20484,17))
            for meg in list_all:
                ii=ii+1;
                print ii
                	
                ### Concrete
                fname_cncrt = data_path + meg + 'Morphed_ico_SemDec_evenodd_ica_equalized_'+event1+'_Mlttap_top10verts_ctf_Coh_'+band+'_150_450_200ms_' + seed_name                   
                stc_cncrt = mne.read_source_estimate(fname_cncrt); #stc_cncrt.data[stc_cncrt.data<0]=0 
                event1_intarget=stc_cncrt.in_label(target_label)
                fname_abs = data_path + meg + 'Morphed_ico_SemDec_evenodd_ica_equalized_'+event2+'_Mlttap_top10verts_ctf_Coh_'+band+'_150_450_200ms_' + seed_name
                stc_abs = mne.read_source_estimate(fname_abs);# stc_abs.data[stc_abs.data<0]=0 
                event2_intarget=stc_abs.in_label(target_label)
                stc_data_subtract=np.subtract(np.mean(np.abs(event1_intarget.data[:,:]),axis=0),np.mean(np.abs(event2_intarget.data[:,:]),axis=0))
#                stc_data_avg=(np.abs(stc_cncrt.data)+np.abs(stc_abs.data))/2.
                X[ii,:,tcnt,bandc]=stc_data_subtract#stc_subtract.data
#                Xa[:,:,ii]=stc_data_avg
            		    		
            for cnt in range(X.shape[1]):
                print cnt
                MT[:,cnt,tcnt,bandc]= mne.stats.ttest_1samp_no_p(X[:,cnt,tcnt,bandc], sigma=1e-3, method='relative')
#            Xavg=np.squeeze(np.mean(X,axis=0))
#            Xav=np.squeeze(np.mean(Xa,axis=2))
            #X = np.transpose(X, [2, 1, 0])	
                
#            MT_stc = mne.SourceEstimate(MT, vertices=stc_abs.vertices, tmin= 0.05, tstep= 0.1, subject='fsaverage')
#            data_out = out_path +  'UVttest_Morphed_ico_SD_ica_evenodd_Sub_Mlttap_'+event1+'_'+event2+'_ctf_Coh_'+band+'_150_450_200ms_' + seed_name_sv+'_'+target_name
#            MT_stc.save(data_out)
                
#            Xavg_stc = mne.SourceEstimate(Xavg, vertices=stc_abs.vertices, tmin= 0.05, tstep= 0.1, subject='fsaverage')
#            data_out = avg_path +  'GrandAverage_Morphed_ico_SD_ica_evenodd_Sub_Mlttap_'+event1+'_'+event2+'_ctf_Coh_'+band+'_150_450_200ms_' + seed_name_sv+'_'+target_name
#            Xavg_stc.save(data_out)
#            Xav_stc = mne.SourceEstimate(Xav, vertices=stc_abs.vertices, tmin= 0.05, tstep= 0.1, subject='fsaverage')
#            data_out = avg_path +  'GrandAverage_Morphed_ico_SD_ica_evenodd_Avg_Mlttap_'+event1+'_'+event2+'_ctf_Coh_'+band+'_150_450_200ms_' + seed_name_sv+'_'+target_name
#            Xav_stc.save(data_out)
#            
pval=stats.t.sf(np.abs(MT),18)*2  
print MT[0,:-1,:,:]   
print pval[0,:-1,:,:]       
#print MT
