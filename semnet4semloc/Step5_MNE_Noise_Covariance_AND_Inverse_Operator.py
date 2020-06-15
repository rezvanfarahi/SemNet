
"""
=========================================================
Noise covariance matrix and inverse operator
=========================================================


"""
# Authors: Alexandre Gramfort <gramfort@nmr.mgh.harvard.edu>
#
# License: BSD (3-clause)

print (__doc__)

import sys
sys.path.insert(1,'/imaging/local/software/mne_python/latest_v0.11')
sys.path.insert(1,'/imaging/rf02/scikit-learn-0.15.0')
sys.path.insert(1,'/home/rf02/.local/lib/python2.7/site-packages')
sys.path.insert(1,'/imaging/local/software/anaconda/2.4.1/2/lib/python2.7/site-packages/scipy')
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
import warnings
warnings.filterwarnings("ignore")
#import matplotlib.pyplot as plt

###############################################################################
data_path = '/imaging/rf02/Semnet/' # root directory for your MEG data
orig_path=data_path
subjects_dir = '/imaging/rf02/Semnet/'    # where your MRI subdirectories are
os.chdir('/imaging/rf02/Semnet/')

# get indices for subjects to be processed from command line input
# 
print (sys.argv)
subject_inds = []#
for ss in sys.argv[1:]:
    subject_inds.append( int( ss ) )

subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]#
print ("subject_inds:")
print (subject_inds)
print ("No rejection")
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


ll = []
for ss in subject_inds:
 ll.append(list_all[ss])

print ("subjects:")
print (ll)



#subjects=[subjects[0]]
for index, meg in enumerate(ll):

   print (meg )#################################################################################
   #read concatenated events and raw file from disc
   
   

   
   print ("reading raw file")
   raw_fname = data_path + meg+ 'block_LD_tsss_filtSL_1_48_ica_out_raw.fif'#'block_LD_tsss_filtSL_1_48_ica_raw.fif'#'SemDec_blocks_tsss_filtSL_1_48_ica_raw.fif'
   raw = mne.io.Raw(raw_fname, preload=True)
   print (raw.info['bads'])
   print ("Reading events"  ) 
   events_fname = orig_path + meg + 'block_LD_tsss_filtSL_1_48_ica_out_raw-eve.fif'#'SemDec_blocks_tsss_filtSL_1_48_ica_raw-eve.fif'
   events = mne.read_events(events_fname)

   print ("Reading evevokeds" )  
   #fname_evoked = directory + 'semloc_ssstf_fft_1_30_clean_ica_raw-ave.fif'
   #evoked = mne.read_evokeds(fname_evoked, condition=0, baseline=(None, 0))
   
#   stim_delay = 0.034 # delay in s
#   events[:,0] = events[:,0]+np.round( raw.info['sfreq']*stim_delay )

   ##################################################################################

   #extracting epochs
   tmin, tmax = -0.3, 0.6
   event_id = { 'Neutral': 4, 'Emotional': 5,'Concrete': 93}
 
   exclude = raw.info['bads']
   include = []
   picks = mne.pick_types(raw.info, meg=True, eeg=True, eog=False, stim=False, include=include, exclude=exclude)
  
   print ("Epoching")
   if subject_inds[index] in np.array([2,3,5]):#np.array([4,8,9]):
        reject = dict(eeg=150e-6, grad=200e-12, mag=4e-12)
   else:
        reject = dict(eeg=120e-6, grad=200e-12, mag=4e-12)#, eog=150e-6)
   epochs_pre = mne.Epochs(raw, events, event_id, tmin, tmax, picks=picks, proj=True, baseline=(None, 0), reject=reject,preload=True)
   print (len(epochs_pre))
   for eegcnt in range(71):
        if eegcnt<10:
            thiseegbad=sum(epochs_pre.drop_log,[]).count('EEG00'+str(eegcnt))
            if thiseegbad>=30:
                raw.info['bads'].append(u'EEG00'+str(eegcnt))                
        else:
            thiseegbad=sum(epochs_pre.drop_log,[]).count('EEG0'+str(eegcnt))
            if thiseegbad>=30:
                raw.info['bads'].append(u'EEG0'+str(eegcnt))
   print (raw.info['bads'])
   picks = mne.pick_types(raw.info, eeg=True, meg=True, eog=True, exclude='bads')
   epochs = mne.Epochs(raw, events, event_id, tmin, tmax, picks=picks,
    		             baseline=(None, 0), proj=True, reject=reject, preload=True)
   print(epochs.drop_log_stats)

   evoked=epochs.average()

   #epochs=mne.epochs.equalize_epoch_counts(epochs)
   ##################################################################################

   #read forward solution
   fname_fwd = data_path + meg + 'ico5_forward_5-3L-EMEG-fwd.fif'

   print ("reading forward solution")
   forward_meeg = mne.read_forward_solution(fname_fwd, surf_ori=True)

   print ("making noise covariance matrix")
   
    
   method_params=dict(diagonal_fixed=dict(grad=0.1, mag=0.1, eeg=0.1))
#   method_params=dict(empirical= dict(store_precision=False, assume_centered= True),diagonal_fixed=dict(grad=0.05, mag=0.05, eeg=0.05,store_precision=False,assume_centered=True),shrunk=dict(shrinkage=np.logspace(-4, 0, 30),store_precision= False, assume_centered= True))

#   noise_cov = mne.compute_covariance(epochs, method=['shrunk','empirical', 'diagonal_fixed'],  tmin=None, tmax=0.,return_estimators=True)#  method_params=method_params)#return_estimators=True,
   noise_cov = mne.compute_covariance(epochs, method=['diagonal_fixed'],  tmin=None, tmax=0., method_params=method_params)#, projs=raw.info['projs'])

   fname_cov = data_path + meg + 'block_LD_tsss_filtSL_1_48_ica_oldreg_noise-cov.fif'
#   noise_cov = mne.read_cov(fname_cov)
   mne.write_cov(fname_cov,noise_cov)
   
   info = evoked.info
   
   print ("making inverse operator")
   inverse_operator_meeg = make_inverse_operator(info, forward_meeg, noise_cov,loose=0.2, depth=None)
   #noise_cov = mne.cov.regularize(noise_cov, epochs.info, mag=0.1, grad=0.1, eeg=0.1, proj=True)

   fname_inv = data_path + meg + 'InvOp_ico5newreg_fft_1_48_clean_ica_EMEG-inv.fif'
   inverse_operator_old = read_inverse_operator(fname_inv)
   sub_src=inverse_operator_old['src']#inverse_operator_meeg['src']
#   mne.add_source_space_distances(sub_src)
   inverse_operator_meeg['src']=sub_src
   #Write the inverse operator
   write_inverse_operator(data_path + meg + 'InvOp_LD_ico5oldreg_fftSL_1_48_clean_ica_EMEG-inv.fif', inverse_operator_meeg)
#   write_inverse_operator(data_path + subjects[index] + 'InverseOperator_fft_1_30_clean_ica_MEG-inv.fif', inverse_operator_meg)
#   write_inverse_operator(data_path + subjects[index] + 'InverseOperator_fft_1_30_clean_ica_EEG-inv.fif', inverse_operator_eeg)
 
#
#   ###end loop across subjects
