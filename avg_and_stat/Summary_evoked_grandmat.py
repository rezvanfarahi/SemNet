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
sys.path.append('/imaging/local/software/mne_python/latest_v0.11')

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
data_path = '/imaging/rf02/Semnet/' # root directory for your MEG data
os.chdir(data_path)
subjects_dir = '/imaging/rf02/Semnet/'    # where your MRI subdirectories are
# where event-files are
#out_path = '/imaging/rf02/Semnet/stc/permutation/power/final_allbands/abs/' #
#if not os.path.exists(out_path):
#    os.makedirs(out_path)

#label_path = '/imaging/rf02/TypLexMEG/fsaverage/label'

inv_path = '/imaging/rf02/Sement/'
#inv_fname = 'InverseOperator_EMEG-inv.fif'

# get indices for subjects to be processed from command line input
# 
print sys.argv
p_inds = [0]
for ss in sys.argv[1:]:
   p_inds.append( int( ss ) )

subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18] # removed
#p_inds=[0,1,2,3,4,5,6,7,8,9,10]
p_list=[0.1]#,0.05,0.045,0.04,0.03,0.025,0.01,0.008,0.005,0.002,0.001,0.0005,0.0001]
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
tmin, tmax = -0.3, 0.6
# frequency bands with variable number of cycles for wavelets
stim_delay = 0.034 # delayfrequencies<=20))] = 3

event_names=['Visual','Hear','Hand','Pwordc']

ntimes=901
nbands=5
nverts=20484
nsubjs=19
nevents=len(event_names)
#X=np.zeros((nverts,ntimes*nbands,nevents,nsubjs))
Xe=np.zeros((nverts,ntimes,nbands,nevents,nsubjs))
Xo=np.zeros((nverts,ntimes,nbands,nevents,nsubjs))
Xa=np.zeros((nverts,ntimes,nbands,nevents,nsubjs))
#X_beta=np.zeros((20484,16))
bands=['delta','theta','alpha','beta','gamma']#'theta','alpha','beta',
for p_threshold in ll:
    ii=-1
    
    for meg in list_all:
        ii=ii+1;
        print ii
        for bcnt,b in enumerate(bands):
            print b
            for eventc, this_event in enumerate(event_names):    
                print this_event
                fname_even = data_path + meg + 'firstMorphed_ico_SemDec_Even_pnt1_48ica_'+this_event+'_Source_Evoked_m300_600_'+b#'mineMorphed_ico_TypLex_icaclean_word_Power_ratio_equalized_m500_700_' + b
                stc_even = mne.read_source_estimate(fname_even)
                Xe[:,:,bcnt,eventc,ii]=stc_even.data
                fname_odd = data_path + meg + 'firstMorphed_ico_SemDec_Odd_pnt1_48ica_'+this_event+'_Source_Evoked_m300_600_'+b#'mineMorphed_ico_TypLex_icaclean_word_Power_ratio_equalized_m500_700_' + b
                stc_odd = mne.read_source_estimate(fname_odd)
                Xo[:,:,bcnt,eventc,ii]=stc_odd.data
#                fname_all = data_path + meg + 'mineMorphed_bandnoisecov_ico_SemDec_ica_'+this_event+'_Power_normori_ratio_equalized_m500_700_200hz_' + b#'mineMorphed_ico_TypLex_icaclean_word_Power_ratio_equalized_m500_700_' + b
#                stc_all = mne.read_source_estimate(fname_all)  
#                Xa[:,:,bcnt,eventc,ii]=np.log(stc_all.data)
                
      ####for PCA!!!
#    for cnt, meg in enumerate(list_all):
#        print cnt
#        path_mtx=data_path+meg+'TF_logratio_even.mat'
#        thisXe=Xe[:,:,:,:,cnt].copy().squeeze()
#        scio.savemat(path_mtx,{'thisXe':thisXe})
#        path_mtx=data_path+meg+'TF_logratio_odd.mat'
#        thisXo=Xo[:,:,:,:,cnt].copy().squeeze()
#        scio.savemat(path_mtx,{'thisXo':thisXo})
#        path_mtx=data_path+meg+'TF_logratio_all.mat'
#        thisXa=Xa[:,:,:,:,cnt].copy().squeeze()
#        scio.savemat(path_mtx,{'thisXa':thisXa})   
#R=np.zeros((nverts,nbands,nevents,nsubjs))
#for ii in range(Xo.shape[0]): #verts
#    for jj in range(Xo.shape[2]):#bands
#        for kk in range(Xo.shape[3]): #events
#            for mm in range(Xo.shape[4]):#subs
#                R[ii,jj,kk,mm]=np.corrcoef(Xe[ii,60:145,jj,kk,mm],Xo[ii,60:145,jj,kk,mm])[0,1]
R=np.zeros((ntimes,nbands,nevents,nsubjs))
for ii in range(Xo.shape[1]): #verts
#    print ii
    for jj in range(Xo.shape[2]):#bands
        for kk in range(Xo.shape[3]): #events
            for mm in range(Xo.shape[4]):#subs
                R[ii,jj,kk,mm]=np.corrcoef(np.exp(Xe[:,ii,jj,kk,mm]),np.exp(Xo[:,ii,jj,kk,mm]))[0,1]     
     
#R3=np.zeros((ntimes,nevents,nsubjs))
#for ii in range(Xo.shape[1]): #times
#        for kk in range(Xo.shape[3]): #events
#            for mm in range(Xo.shape[4]):#subs
#                R3[ii,kk,mm]=np.corrcoef(Xe[:,ii,0,kk,mm],Xo[:,ii,0,kk,mm])[0,1]
#     
#     