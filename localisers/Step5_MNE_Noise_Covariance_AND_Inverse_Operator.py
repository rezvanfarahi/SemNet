
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

subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13]
print "subject_inds:"
print subject_inds
print "No rejection"
list_all =  ['/meg16_0045/160303/', #8
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

print "subjects:"
print ll



#subjects=[subjects[0]]
for index, meg in enumerate(ll):

   print meg #################################################################################
   #read concatenated events and raw file from disc
   
   

   
   print "reading raw file"
   raw_fname = data_path + meg + 'block_localisers_tsss_filt_pnt1_48_ica_raw.fif'#'block_milk_tsss_filt_pnt1_48_ica_raw.fif' #clean_ssp_
   raw = mne.io.Raw(raw_fname, preload=True)

   print "Reading events"   
   events_fname = data_path + meg + 'block_localisers_tsss_ds_raw-eve.fif'
   events = mne.read_events(events_fname)

   print "Reading evevokeds"   
   #fname_evoked = directory + 'semloc_ssstf_fft_1_48_clean_ica_raw-ave.fif'
   #evoked = mne.read_evokeds(fname_evoked, condition=0, baseline=(None, 0))
   
   stim_delay = 0.034 # delay in s
   events[:,0] = events[:,0]+np.round( raw.info['sfreq']*stim_delay )

   ##################################################################################

   #extracting epochs
   print "epochs"
   tmin=-0.1
   tmax=0.6
   picks = mne.pick_types(raw.info, eeg=True, meg=True, eog=True, exclude='bads')
   event_id = {'audio': 1,  'colour': 3, 'grey': 4, 'shapes': 7, 'scrambled': 8}#, 'button': 4096}
   color = {1: 'blue', 3: 'red', 4: 'green', 7: 'c', 8: 'black'}#, 4096: 'yellow'}
   if subject_inds[index] in np.array([0]):#np.array([4,8,9]):
        reject = dict(eeg=150e-6, grad=200e-12, mag=4e-12)#, eog=150e-6)
   else:
        reject = dict(eeg=120e-6, grad=200e-12, mag=4e-12)#, eog=150e-6)
    
   epochs = mne.Epochs(raw, events, event_id, tmin, tmax, picks=picks,
    		             baseline=(None,0), proj=True, reject=reject, preload=True)
   
   evoked=epochs.average()
   #epochs=mne.epochs.equalize_epoch_counts(epochs)
   ##################################################################################

   #read forward solution
   fname_fwd = data_path + meg + 'ico5_forward_5-3L-EMEG-fwd.fif'

   print "reading forward solution"
   forward_meeg = mne.read_forward_solution(fname_fwd, surf_ori=True)
#   forward_meg = mne.pick_types_forward(forward_meeg, meg=True, eeg=False)
   # Alternatively, you can just load a forward solution that is restricted
#   forward_eeg = mne.pick_types_forward(forward_meeg, meg=False, eeg=True)
   print "making noise covariance matrix"
   


   #Compute the noise covariance on baseline
    
#   method_params=dict(diagonal_fixed=dict(grad=0.1, mag=0.1, eeg=0.1))
#   method_params=dict(empirical= dict(store_precision=False, assume_centered= True),diagonal_fixed=dict(grad=0.05, mag=0.05, eeg=0.05,store_precision=False,assume_centered=True),shrunk=dict(shrinkage=np.logspace(-4, 0, 30),store_precision= False, assume_centered= True))

   noise_cov = mne.compute_covariance(epochs, method=['shrunk','empirical', 'diagonal_fixed'],  tmin=None, tmax=0.)#  method_params=method_params)#return_estimators=True,
#   noise_cov = mne.compute_covariance(epochs, method=['diagonal_fixed'],  tmin=None, tmax=0., method_params=method_params, projs=raw.info['projs'])

   fname_cov = data_path + meg + 'localisers_tsss_filt_pnt1_48_ica_newreg_noise-cov.fif'#'localisers_tsss_filt_pnt1_48_ica_oldreg_noise-cov.fif'
   noise_cov = mne.read_cov(fname_cov)
#   mne.write_cov(fname_cov,noise_cov)
#   noise_cov_meg = mne.compute_covariance(epochs_meg, method=['empirical','shrunk', 'diagonal_fixed'], tmin=None, tmax=0.)
#   noise_cov_eeg = mne.compute_covariance(epochs_eeg, method=['empirical','shrunk', 'diagonal_fixed'], tmin=None, tmax=0.)
#   regularize noise cov...
#   ####noise_cov = mne.cov.regularize(noise_cov, epochs.info, mag=0.1, grad=0.1, eeg=0.1, proj=True) ##I checked with and without proj, both exactly the same in this case
#   mne.viz.plot_cov(noise_cov, raw.info, colorbar=True, proj=True)
#   make an M/EEG, MEG-only, and EEG-only inverse operators
   ###noise_cov_eeg = mne.cov.regularize(noise_cov_eeg, epochs_eeg.info, eeg=0.1, proj=True)
   ###noise_cov_meg = mne.cov.regularize(noise_cov_meg, epochs_meg.info, mag=0.1, grad=0.1, proj=True)
   
   
   info = evoked.info
   #noise_cov = mne.compute_covariance(epochs, tmin=None, tmax=0.)
#   info_eeg = evoked_eeg.info
#   info_meg = evoked_meg.info
   print "a"
   inverse_operator_meeg = make_inverse_operator(info, forward_meeg, noise_cov,loose=0.2, depth=None)
   #noise_cov = mne.cov.regularize(noise_cov, epochs.info, mag=0.1, grad=0.1, eeg=0.1, proj=True)
#   print "b"
#   inverse_operator_meg = make_inverse_operator(info_meg, forward_meg, noise_cov_meg, loose=0.2, depth=None)
#   print "c"
#   inverse_operator_eeg = make_invrse_operator(info_eeg, forward_eeg, noise_cov_eeg, loose=0.2, depth=None)
   fname_inv = data_path + meg + 'InvOp_ico5oldreg_fft_pnt1_48_localisers_ica_EMEG-inv.fif'
   inverse_operator_old = read_inverse_operator(fname_inv)
   sub_src=inverse_operator_old['src']#nverse_operator_meeg['src']
#   mne.add_source_space_distances(sub_src)
   inverse_operator_meeg['src']=sub_src
   #Write the inverse operator
   write_inverse_operator(data_path + meg + 'InvOp_ico5newreg_fft_pnt1_48_localisers_ica_EMEG-inv.fif', inverse_operator_meeg)
#   write_inverse_operator(data_path + subjects[index] + 'InverseOperator_fft_1_48_clean_ica_MEG-inv.fif', inverse_operator_meg)
#   write_inverse_operator(data_path + subjects[index] + 'InverseOperator_fft_1_48_clean_ica_EEG-inv.fif', inverse_operator_eeg)
 

   ###end loop across subjects