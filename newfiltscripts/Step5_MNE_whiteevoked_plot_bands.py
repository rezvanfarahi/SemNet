# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 02:43:01 2017

@author: rf02
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 01:59:03 2017

@author: rf02
"""


"""
=========================================================
Noise covariance matrix and inverse operator
=========================================================


"""
# Authors: Alexandre Gramfort <gramfort@nmr.mgh.harvard.edu>
#
# License: BSD (3-clause)

print __doc__

import sys
sys.path.insert(1,'/imaging/local/software/mne_python/latest_v0.15')
sys.path.insert(1,'/imaging/rf02/scikit-learn-0.15.0')
sys.path.insert(1,'/home/rf02/.local/lib/python2.7/site-packages')
#sys.path.insert(1,'/imaging/local/software/anaconda/2.4.1/2/lib/python2.7/site-packages/scipy')
#import scipy
#print scipy.__version__
# for qsub
# (add ! here if needed) /imaging/local/software/anaconda/latest/x86_64/bin/python
#sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
###

import numpy as np
import mne

from mne.minimum_norm import ( make_inverse_operator, write_inverse_operator, read_inverse_operator)
import os
import matplotlib.pyplot as plt
from mne.io.pick import channel_indices_by_type

###############################################################################
data_path = '/imaging/rf02/Semnet/' # root directory for your MEG data
orig_path=data_path
subjects_dir = '/imaging/rf02/Semnet/'    # where your MRI subdirectories are
os.chdir('/imaging/rf02/Semnet/')

# get indices for subjects to be processed from command line input
# 
print sys.argv
subject_inds = []#
for ss in sys.argv[1:]:
    subject_inds.append( int( ss ) )

subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]#
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


ll = []
for ss in subject_inds:
 ll.append(list_all[ss])

print "subjects:"
print ll


bands=['theta','alpha','beta','gamma']
freqsmin=np.array([4,8,12,30])#np.array([4,8,30])
freqsmax=np.array([8,12,30,45])#np.array([8,12,48])
#subjects=[subjects[0]]
for index, meg in enumerate(ll):

   print meg; subject=subject_inds[index]; #################################################################################
   #read concatenated events and raw file from disc
   
   

#   
#   print "reading raw file"
#   raw_fname = data_path + meg+ 'SemDec_blocks_tsss_filtnew_pnt1_48_ica_raw.fif'
#   raw = mne.io.Raw(raw_fname, preload=True)
#   print raw.info['bads']
#   print "Reading events"   
#   events_fname = orig_path + meg + 'SemDec_blocks_tsss_filtnew_pnt1_48_ica_raw-eve.fif'
#   events = mne.read_events(events_fname)
#
#   print "Reading evevokeds"   
#   #fname_evoked = directory + 'semloc_ssstf_fft_1_30_clean_ica_raw-ave.fif'
#   #evoked = mne.read_evokeds(fname_evoked, condition=0, baseline=(None, 0))
#   
##   stim_delay = 0.034 # delay in s
##   events[:,0] = events[:,0]+np.round( raw.info['sfreq']*stim_delay )
#
#   ##################################################################################
#
#   #extracting epochs for the original unfiltered signal
#   tmin, tmax = -0.3, 0.6
#   event_id = {'visual': 1, 'hear': 2, 'hand': 3, 'pwordc': 6}#'neutral': 4, 'emotional': 5,
# 
#   exclude = raw.info['bads']
#   include = []
#   picks = mne.pick_types(raw.info, meg=True, eeg=True, eog=False, stim=False, include=include, exclude=exclude)
#   print "Epoching"
#   if subject in np.array([2,3,5,18]):#np.array([4,8,9]):
#        reject = dict(eeg=150e-6, grad=200e-12, mag=4e-12)
#   else:
#        reject = dict(eeg=120e-6, grad=200e-12, mag=4e-12)#, eog=150e-6)
#   epochspre = mne.Epochs(raw, events, event_id, tmin, tmax, picks=picks, proj=True, baseline=(None, 0), reject=None,preload=True)
#   chidx=channel_indices_by_type(epochspre.info)
#   
#   drop_inds=list()
#   edata=epochspre.get_data()
#   for epii in range(edata.shape[0]):
#       for ch_type in ['grad','mag','eeg']:
#           for chjj in chidx[ch_type]:
#               thresh=reject[ch_type]
#               if np.max(edata[epii,chjj,:])-np.min(edata[epii,chjj,:])>thresh:
#                   drop_inds.append(epii)
#   drop_inds=np.unique(np.asarray(drop_inds))
#   print len(epochspre)/600.
#   raw_orig=raw.copy()
   for bcnt,b in enumerate(bands):
       print b
#       raw=raw_orig.copy()
#       freqmin=freqsmin[bcnt].copy();freqmax=freqsmax[bcnt].copy()
##       raw.filter(l_freq=freqmin, h_freq=freqmax,  l_trans_bandwidth=0.05,h_trans_bandwidth=0.1, picks=picks, method='fft',filter_length='50s')
#       if b in ['theta','alpha','beta']:
#           raw.filter(l_freq=freqmin, h_freq=freqmax,  l_trans_bandwidth='auto',h_trans_bandwidth='auto', picks=picks, method='fir',filter_length='auto')#,phase='zero-double')
#       else:
#           raw.filter(l_freq=freqmin,  h_freq=None, l_trans_bandwidth=5, picks=picks, method='fir',filter_length='auto')#,phase='zero-double')
#
#       epochs = mne.Epochs(raw, events, event_id, tmin, tmax, picks=picks, proj=True, baseline=(None, 0), reject=None,preload=True)
#       epochs.drop(drop_inds)
#       print len(epochs)/600.
#       
#    #   for eegcnt in range(71):
#    #        if eegcnt<10:
#    #            thiseegbad=sum(epochs.drop_log,[]).count('EEG00'+str(eegcnt))
#    #            if thiseegbad>=90:
#    #                raw.info['bads'].append(u'EEG00'+str(eegcnt))                
#    #        else:
#    #            thiseegbad=sum(epochs.drop_log,[]).count('EEG0'+str(eegcnt))
#    #            if thiseegbad>=90:
#    #                raw.info['bads'].append(u'EEG0'+str(eegcnt))
#       evoked=epochs.average()
#    
#       #epochs=mne.epochs.equalize_epoch_counts(epochs)
#       ##################################################################################
#    
#       #read forward solution
#       fname_fwd = data_path + meg + 'ico5_forward_5-3L-EMEG-fwd.fif'
#    
#       print "reading forward solution"
#       forward_meeg = mne.read_forward_solution(fname_fwd, surf_ori=True)
#    #   forward_meg = mne.pick_types_forward(forward_meeg, meg=True, eeg=False)
#       # Alternatively, you can just load a forward solution that is restricted
#    #   forward_eeg = mne.pick_types_forward(forward_meeg, meg=False, eeg=True)
#       print "making noise covariance matrix"
       
    
    
       #Compute the noise covariance on baseline
        
    #   method_params=dict(diagonal_fixed=dict(grad=0.1, mag=0.1, eeg=0.1))
    #   method_params=dict(empirical= dict(store_precision=False, assume_centered= True),diagonal_fixed=dict(grad=0.05, mag=0.05, eeg=0.05,store_precision=False,assume_centered=True),shrunk=dict(shrinkage=np.logspace(-4, 0, 30),store_precision= False, assume_centered= True))
    
#       noise_cov = mne.compute_covariance(epochs, method=['shrunk','empirical', 'diagonal_fixed'],  tmin=None, tmax=0.,return_estimators=True)#  method_params=method_params)#return_estimators=True,
    #   noise_cov = mne.compute_covariance(epochs, method=['diagonal_fixed'],  tmin=None, tmax=0., method_params=method_params, projs=raw.info['projs'])
    
       fname_cov = data_path + meg + 'SemDec_blocks_tsss_filtnew_pnt1_48_ica_newreg_'+b+'_noise-cov.fif'
       noise_cov = mne.read_cov(fname_cov)
#       mne.write_cov(fname_cov,noise_cov[0])
       figname_in='/imaging/rf02/Semnet/results/evoked/'+b+'_noisecov_filtnew/subj'+str(subject)+'_evoked-ave.fif'#.jpg
       evoked=mne.read_evokeds(fname=figname_in,baseline=(None,0),kind='average',proj=True,allow_maxshield=True)
       fig=evoked[0].plot(show=False)
       figname_out='/imaging/rf02/Semnet/results/evoked/'+b+'_noisecov_filtnew/subj'+str(subject)+'.jpg'

#       evoked.save(figname_out)
       fig.savefig(figname_out)
       plt.close("all")
##   noise_cov_meg = mne.compute_covariance(epochs_meg, method=['empirical','shrunk', 'diagonal_fixed'], tmin=None, tmax=0.)
##   noise_cov_eeg = mne.compute_covariance(epochs_eeg, method=['empirical','shrunk', 'diagonal_fixed'], tmin=None, tmax=0.)
##   regularize noise cov...
##   ####noise_cov = mne.cov.regularize(noise_cov, epochs.info, mag=0.1, grad=0.1, eeg=0.1, proj=True) ##I checked with and without proj, both exactly the same in this case
##   mne.viz.plot_cov(noise_cov, raw.info, colorbar=True, proj=True)
##   make an M/EEG, MEG-only, and EEG-only inverse operators
#   ###noise_cov_eeg = mne.cov.regularize(noise_cov_eeg, epochs_eeg.info, eeg=0.1, proj=True)
#   ###noise_cov_meg = mne.cov.regularize(noise_cov_meg, epochs_meg.info, mag=0.1, grad=0.1, proj=True)
#   
#   
#       info = evoked.info
#       
#       #noise_cov = mne.compute_covariance(epochs, tmin=None, tmax=0.)
#    #   info_eeg = evoked_eeg.info
#    #   info_meg = evoked_meg.info
#       print "a"
#       inverse_operator_meeg = make_inverse_operator(info, forward_meeg, noise_cov,loose=0.2, depth=None)
#       #noise_cov = mne.cov.regularize(noise_cov, epochs.info, mag=0.1, grad=0.1, eeg=0.1, proj=True)
#    #   print "b"
#    #   inverse_operator_meg = make_inverse_operator(info_meg, forward_meg, noise_cov_meg, loose=0.2, depth=None)
#    #   print "c"
#    #   inverse_operator_eeg = make_invrse_operator(info_eeg, forward_eeg, noise_cov_eeg, loose=0.2, depth=None)
#       if subject in np.array([18]):
#           fname_inv = data_path + meg + 'InvOp_ico5newreg_filtnew_pnt1_30_clean_ica_EMEG-inv.fif'
#           inverse_operator_old = read_inverse_operator(fname_inv)
#           sub_src=inverse_operator_old['src']#inverse_operator_meeg['src']
#       else:
#           fname_inv = data_path + meg + 'InvOp_ico5newreg_fft_1_48_clean_ica_EMEG-inv.fif'
#           inverse_operator_old = read_inverse_operator(fname_inv)
#           sub_src=inverse_operator_old['src']
#        #   mne.add_source_space_distances(sub_src)
#       inverse_operator_meeg['src']=sub_src
#       #Write the inverse operator
#       write_inverse_operator(data_path + meg + 'InvOp2_ico5newreg_filtnew_pnt1_48_clean_ica_'+b+'_EMEG-inv.fif', inverse_operator_meeg)
##       write_inverse_operator(data_path + subjects[index] + 'InverseOperator_fft_1_48_clean_ica_MEG-inv.fif', inverse_operator_meg)
##       write_inverse_operator(data_path + subjects[index] + 'InverseOperator_fft_1_48_clean_ica_EEG-inv.fif', inverse_operator_eeg)
#     
#
#   ##end loop across subjects
