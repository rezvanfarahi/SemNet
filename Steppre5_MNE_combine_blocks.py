
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
sys.path.insert(1,'/imaging/local/software/mne_python/latest_v0.9full')
sys.path.insert(1,'/imaging/rf02/scikit-learn-0.15.0')
sys.path.insert(1,'/home/rf02/.local/lib/python2.7/site-packages')

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
orig_path='/group/erp/data/olaf.hauk/MEG/TypLexM/data/MF22_new/'
subjects_dir = '/imaging/rf02/Semnet/'    # where your MRI subdirectories are
os.chdir('/imaging/rf02/Semnet/')

# get indices for subjects to be processed from command line input
# 
print sys.argv
subject_inds = []#
for ss in sys.argv[1:]:
    subject_inds.append( int( ss ) )

subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]
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
            '/meg16_0123/160708/', #23 
            '/meg16_0125/160712/', #24 
            ]


ll = []
for ss in subject_inds:
 ll.append(list_all[ss])

print "subjects:"
print ll



#subjects=[subjects[0]]
for index, meg in enumerate(ll):

   print index #################################################################################
   #read concatenated events and raw file from disc
      

   
   print "reading raw files"
   raw_fname1 = data_path+meg+'block_fruit_tsss_filtnew_pnt1_30_ica_raw.fif'#directory + 'semloc_ssstf_fft_1_30_clean_ica_raw.fif'
   raw1 = mne.io.Raw(raw_fname1, preload=True)
   raw_fname2 = data_path+meg+'block_milk_tsss_filtnew_pnt1_30_ica_raw.fif'#directory + 'semloc_ssstf_fft_1_30_clean_ica_raw.fif'
   raw2 = mne.io.Raw(raw_fname2, preload=True)
   raw_fname3 = data_path+meg+'block_odour_tsss_filtnew_pnt1_30_ica_raw.fif'#directory + 'semloc_ssstf_fft_1_30_clean_ica_raw.fif'
   raw3 = mne.io.Raw(raw_fname3, preload=True)

   print "Reading events"   
   stim_delay = 0.034 # delay in s
   events_fname1 = data_path + meg + 'block_fruit_tsss_raw-eve.fif'#orig_path + meg + 'semloc_raw_ssstnew.txt'
   events1 = mne.read_events(events_fname1)
   events1[:,0]= np.round( raw1.info['sfreq']*stim_delay )+events1[:,0]
   events_fname2 = data_path + meg + 'block_milk_tsss_raw-eve.fif'#orig_path + meg + 'semloc_raw_ssstnew.txt'
   events2 = mne.read_events(events_fname2)
   events2[:,0]= np.round( raw2.info['sfreq']*stim_delay )+events2[:,0]
   events_fname3 = data_path + meg + 'block_odour_tsss_raw-eve.fif'#orig_path + meg + 'semloc_raw_ssstnew.txt'
   events3 = mne.read_events(events_fname3)
   events3[:,0]= np.round( raw3.info['sfreq']*stim_delay )+events3[:,0]


   print "Reading evevokeds"   
   #fname_evoked = directory + 'semloc_ssstf_fft_1_30_clean_ica_raw-ave.fif'
   #evoked = mne.read_evokeds(fname_evoked, condition=0, baseline=(None, 0))
   raws=[raw1,raw2,raw3]
   
   eventss=[events1,events2,events3]
   raw_fname = data_path+meg+'SemDec_blocks_tsss_filtnew_pnt1_30_ica_raw.fif'
   raw,events=mne.concatenate_raws(raws,events_list=eventss)
   eventsch=events.copy()
    ###############################################################################
    # Read epochs for all channels, removing a bad one
#    raw.plot(scalings=dict(eeg=100e-6))
   print "epochs"
   tmin = -0.3
   tmax = 0.6 
   picks = mne.pick_types(raw.info, eeg=True, meg=True, eog=True, exclude='bads')
   event_id = {'visual': 1, 'hear': 2, 'hand': 3, 'neutral': 4, 'emotional': 5,'pword': 6}
   if subject_inds[index] in np.array([3,5]):#np.array([4,8,9]):
        reject = dict(eeg=120e-6, grad=200e-12, mag=4e-12)
   else:
        reject = dict(eeg=100e-6, grad=200e-12, mag=4e-12)#, eog=150e-6)
    
   epochs = mne.Epochs(raw, eventsch, event_id, tmin, tmax, picks=picks,
    		             baseline=(None, 0), proj=True, reject=reject, preload=True)
   for eegcnt in range(71):
        if eegcnt<10:
            thiseegbad=sum(epochs.drop_log,[]).count('EEG00'+str(eegcnt))
            if thiseegbad>=90:# and eegcnt not in range(8):
                raw.info['bads'].append(u'EEG00'+str(eegcnt))                
        else:
            thiseegbad=sum(epochs.drop_log,[]).count('EEG0'+str(eegcnt))
            if thiseegbad>=90:
                raw.info['bads'].append(u'EEG0'+str(eegcnt))
   print raw.info['bads']
   picks = mne.pick_types(raw.info, eeg=True, meg=True, eog=True, exclude='bads')
   epochs = mne.Epochs(raw, eventsch, event_id, tmin, tmax, picks=picks,
    		             baseline=(None, 0), proj=True, reject=reject, preload=True)
   print len(epochs)/900.
   raw.save(raw_fname, overwrite=True)
   mne.write_events(raw_fname[:-4]+'-eve.fif', events)
   
