
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

from mne import (io, spatial_tris_connectivity, compute_morph_matrix,
                 grade_to_tris)
from mne.stats import (spatio_temporal_cluster_1samp_test,
                       summarize_clusters_stc)


###############################################################################
data_path =  '/imaging/rf02/Semnet/' # root directory for your MEG data
os.chdir(data_path)
subjects_dir =  '/imaging/rf02/Semnet/'    # where your MRI subdirectories are
# where event-files are
out_path = '/imaging/rf02/Semnet/stc/uvttest/connectivity/evenodd/top10verts/'    # where event files are
if not os.path.exists(out_path):
    os.makedirs(out_path)

avg_path = '/imaging/rf02/Semnet/stc/GrandAverage/connectivity/evenodd/top10verts/'    # where event files are
if not os.path.exists(avg_path):
    os.makedirs(avg_path)

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



#'lh.lateralorbitofrontal','rh.lateralorbitofrontal']#,
#labellist=['lhATL_latmed-lh','rhATL_latmed-rh','lhAG_hsmetafinal5-lh','lhposterior_inferiorfrontal-lh','lhmiddle_temporal-lh']#labellist = ['lh.ATL_latmed','lh.supramarginal','lh.posterior_inferiorfrontal','lh.middle_temporal']#labellist = [ 'lh.posterior_fusiform','rh.posterior_fusiform', 'lh.lateraloccipital','rh.lateraloccipital','lh.precentral','rh.precentral', 'lh.ATLsmall','rh.ATLsmall']#labellist = [ 'leftATL-lh','rightATL-rh']#, 'medialtempright-rh','medialtempleft-lh']'lh.postcentral','rh.postcentral',
#labellist_path = [label_path+'0_ATL_latmed_ctf-lh.label',label_path+'1_AG_hsmetafinal5_ctf-lh.label',label_path+'2_col_shape_c5_ctf-lh.label',label_path+'3_col_shape_c5_ctf-rh.label',label_path+'4_auditory_c_ctf-lh.label',label_path+'5_auditory_c_ctf-rh.label',label_path+'6_hand_c_ctf-lh.label',label_path+'7_hand_c_ctf-rh.label']#[ 'lh.posterior_fusiform','rh.posterior_fusiform', 'lh.lateraloccipital','rh.lateraloccipital','lh.precentral','rh.precentral', 'lh.ATL_latmed','rh.ATL_latmed']#labellist = [ 'leftATL-lh','rightATL-rh']#, 'medialtempright-rh','medialtempleft-lh']'lh.postcentral','rh.postcentral',

#labellist1=['lh0_ATL_latmed_ctf', 'lh1_AG_hsmetafinal5_ctf', 'lh2_col_shape_c5_ctf', 'rh3_col_shape_c5_ctf', 'lh4_auditory_c_ctf', 'rh5_auditory_c_ctf', 'lh6_hand_c_ctf', 'rh7_hand_c_ctf']#'lh6_hand_c', 'rh7_hand_c']#'lh4_auditory_c', 'rh5_auditory_c']#, 'hand_c-lh', 'hand_c-rh']
labellist1 = ['lh0_ATL_ventromedial_top10verts_ctf','lh1_ATL_dorsolateral_top10verts_ctf','lh2_AG_anterior_top10verts_ctf','lh3_AG_posterior_top10verts_ctf','lh4_hand_postcentral_top10verts_ctf','lh5_hand_precentral_top10verts_ctf','lh6_col_shape_c5_top10verts_ctf']#,label_path+'5_auditory_c_ctf-rh.label',label_path+'6_hand_c_ctf-lh.label',label_path+'7_hand_c_ctf-rh.label']#,label_path+'2_col_shape_c5_ctf-lh.label',label_path+'3_col_shape_c5_ctf-rh.label',label_path+'4_auditory_c_ctf-lh.label',label_path+'5_auditory_c_ctf-rh.label',label_path+'6_hand_c_ctf-lh.label',label_path+'7_hand_c_ctf-rh.label']#[ 'lh.posterior_fusiform','rh.posterior_fusiform', 'lh.lateraloccipital','rh.lateraloccipital','lh.precentral','rh.precentral', 'lh.ATL_latmed','rh.ATL_latmed']#labellist = [ 'leftATL-lh','rightATL-rh']#, 'medialtempright-rh','medialtempleft-lh']'lh.postcentral','rh.postcentral',

labellist_sv1=['lhATL_ventromedial-lh','lhATL_dorsolateral-lh','lhAG_anterior-lh','lhAG_posterior-lh','lhhnd_postc-lh','lhhnd_prec-lh','lhvis-lh','rhaud-rh','lhhnd-lh','rhhnd-rh']

#labellist = [ 'lh.posterior_fusiform','rh.posterior_fusiform']#,'lh.precentral','rh.precentral', 'lh.lateraloccipital','rh.lateraloccipital']#labellist = [ 'posterior_fusiform-lh','rh.lateralorbitofron-rh', 'rh.medialorbitofron-rh']#labellist = ['leftATL-rh','rightATL-lh']#, 'rightATL-rh']#'atlleft-lh']#, 'atlright-rh']#,'medialtempright-rh','medialtempleft-lh']
conn_methods=['Coherence']#,'DWPLI']
bands=['Theta','Alpha','Beta','Gamma']
labellist=[labellist1[ii] for ii in np.array([0,1,2,3,4,5,6])]
labellist_sv=[labellist_sv1[ii] for ii in np.array([0,1,2,3,4,5,6])]
event1='Visual'
event2='Hear'
ntimes=4
for label_name, label_name_sv in zip(labellist, labellist_sv):
    print label_name
    for con_method in conn_methods:
        for band in bands:
            MT=np.zeros((20484,ntimes))
            ii=-1; print band
            X=np.zeros((19,20484,ntimes))
            Xa=np.zeros((20484,ntimes,19))
            #    	X_theta=np.zeros((20484,17))
            #    	X_alpha=np.zeros((20484,17))
            #    	X_beta=np.zeros((20484,17))
            #    	X_gamma=np.zeros((20484,17))
            for meg in list_all:
                ii=ii+1;
                print ii
                	
                ### Concrete
                fname_cncrt = data_path + meg + 'Morphed_ico_SemDec_evenodd_ica_equalized_'+event1+'_Mlttap_top10verts_ctf_Coh_'+band+'_150_450_200ms_' + label_name                   
                stc_cncrt = mne.read_source_estimate(fname_cncrt); #stc_cncrt.data[stc_cncrt.data<0]=0 
                fname_abs = data_path + meg + 'Morphed_ico_SemDec_evenodd_ica_equalized_'+event2+'_Mlttap_top10verts_ctf_Coh_'+band+'_150_450_200ms_' + label_name
                stc_abs = mne.read_source_estimate(fname_abs);# stc_abs.data[stc_abs.data<0]=0 	
                stc_data_subtract=np.subtract(np.abs(stc_cncrt.data[:,:ntimes]),np.abs(stc_abs.data[:,:ntimes]))
                stc_data_avg=(np.abs(stc_cncrt.data[:,:ntimes])+np.abs(stc_abs.data[:,:ntimes]))/2.
                X[ii,:,:]=stc_data_subtract#stc_subtract.data
                Xa[:,:,ii]=stc_data_avg
            		    		
            for cnt in range(X.shape[2]):
                print cnt
                MT[:,cnt]= mne.stats.ttest_1samp_no_p(X[:,:,cnt], sigma=1e-3, method='relative')
            Xavg=np.squeeze(np.mean(X,axis=0))
            Xav=np.squeeze(np.mean(Xa,axis=2))
            #X = np.transpose(X, [2, 1, 0])	
            MT_stc = mne.SourceEstimate(MT, vertices=stc_abs.vertices, tmin= 0.05, tstep= 0.1, subject='fsaverage')
            data_out = out_path +  'UVttest_Morphed_ico_SD_ica_evenodd_Sub_Mlttap_'+event1+'_'+event2+'_ctf_'+con_method+'_'+band+'_50_550_200ms_lh' + label_name_sv[0:-3] 
            MT_stc.save(data_out)
            Xavg_stc = mne.SourceEstimate(Xavg, vertices=stc_abs.vertices, tmin= 0.05, tstep= 0.1, subject='fsaverage')
            data_out = avg_path +  'GrandAverage_Morphed_ico_SD_ica_evenodd_Sub_Mlttap_'+event1+'_'+event2+'_ctf_'+con_method+'_'+band+'_50_550_200ms_lh' + label_name_sv[0:-3] 
            Xavg_stc.save(data_out)
            Xav_stc = mne.SourceEstimate(Xav, vertices=stc_abs.vertices, tmin= 0.05, tstep= 0.1, subject='fsaverage')
            data_out = avg_path +  'GrandAverage_Morphed_ico_SD_ica_evenodd_Avg_Mlttap_'+event1+'_'+event2+'_ctf_'+con_method+'_'+band+'_50_550_200ms_lh' + label_name_sv[0:-3] 
            Xav_stc.save(data_out)
            
            

