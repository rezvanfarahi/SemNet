
"""
=========================================================
TF representation of typlex for frequency bands in source labels
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
sys.path.append('/imaging/local/software/mne_python/latest_v0.9')
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
data_path = '/imaging/rf02/TypLexMEG/' # root directory for your MEG data
os.chdir(data_path)
subjects_dir = '/home/rf02/rezvan/TypLexMEG/'    # where your MRI subdirectories are
# where event-files are
out_path = '/imaging/rf02/TypLexMEG/icaanalysis_results/typlex/stc/UVttest/connectivity/WB_hubs/newp/'    # where event files are


inv_path = '/imaging/rf02/TypLexMEG/'
inv_fname = 'InverseOperator_EMEG-inv.fif'

# get indices for subjects to be processed from command line input
# 

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


#'lh.lateralorbitofrontal','rh.lateralorbitofrontal']#,
labellist=['ATL_latmed-lh','supramarginal-lh','posterior_inferiorfrontal-lh','middle_temporal-lh']
#labellist = ['lh.ATL_latmed','lh.supramarginal','lh.posterior_inferiorfrontal','lh.middle_temporal']#[ 'lh.posterior_fusiform','rh.posterior_fusiform', 'lh.lateraloccipital','rh.lateraloccipital','lh.precentral','rh.precentral', 'lh.ATL_latmed','rh.ATL_latmed']#labellist = [ 'leftATL-lh','rightATL-rh']#, 'medialtempright-rh','medialtempleft-lh']'lh.postcentral','rh.postcentral',
#labellist = [ 'lh.posterior_fusiform','rh.posterior_fusiform', 'lh.lateraloccipital','rh.lateraloccipital','lh.precentral','rh.precentral', 'lh.ATL_latmed','rh.ATL_latmed']#labellist = [ 'leftATL-lh','rightATL-rh']#, 'medialtempright-rh','medialtempleft-lh']'lh.postcentral','rh.postcentral',
#labellist = [ 'lh.posterior_fusiform','rh.posterior_fusiform']#,'lh.precentral','rh.precentral', 'lh.lateraloccipital','rh.lateraloccipital']#labellist = [ 'posterior_fusiform-lh','rh.lateralorbitofron-rh', 'rh.medialorbitofron-rh']#labellist = ['leftATL-rh','rightATL-lh']#, 'rightATL-rh']#'atlleft-lh']#, 'atlright-rh']#,'medialtempright-rh','medialtempleft-lh']
conn_methods=['WPLI2d','Coherence']#,'DWPLI']
bands=['Theta','Alpha','Beta','Gamma']
for label_name in labellist:
    print label_name
    for con_method in conn_methods:
        for band in bands:
            MT=np.zeros((20484,2))
            ii=-1; print band
            X=np.zeros((17,20484,2))
            #    	X_theta=np.zeros((20484,17))
            #    	X_alpha=np.zeros((20484,17))
            #    	X_beta=np.zeros((20484,17))
            #    	X_gamma=np.zeros((20484,17))
            for meg in list_all:
                ii=ii+1;
                print ii
                	
                ### word
                fname_cncrt = data_path + meg + 'Morphed_ico_typlex_ica_equalized_word_Mlttap_bestV_ctf_'+con_method+'_'+band+'_150_450_200ms_' + label_name[0:-3] 
                  
                stc_cncrt = mne.read_source_estimate(fname_cncrt); stc_cncrt.data[stc_cncrt.data<0]=0 
                fname_abs = data_path + meg + 'Morphed_ico_typlex_ica_equalized_nonword_Mlttap_bestV_ctf_'+con_method+'_'+band+'_150_450_200ms_' + label_name[0:-3] 
                stc_abs = mne.read_source_estimate(fname_abs); stc_abs.data[stc_abs.data<0]=0 	
                stc_data_subtract=np.subtract(np.abs(stc_cncrt.data),np.abs(stc_abs.data))
                X[ii,:,:]=stc_data_subtract#stc_subtract.data
            		    		
            for cnt in range(X.shape[2]):
                print cnt
                MT[:,cnt]= mne.stats.ttest_1samp_no_p(X[:,:,cnt], sigma=1e-3, method='relative')
              
            #X = np.transpose(X, [2, 1, 0])	
            MT_stc = mne.SourceEstimate(MT, vertices=stc_abs.vertices, tmin= 0.15, tstep= 0.1, subject='fsaverage')
            data_out = out_path +  'UVttest_Morphed_ico_typlex_ica_equalized_subtract_Mlttap_bestV_ctf_'+con_method+'_'+band+'_150_450_200ms_lh' + label_name[0:-3] 
            MT_stc.save(data_out)
            
            

