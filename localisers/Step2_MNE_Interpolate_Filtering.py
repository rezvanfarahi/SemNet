"""
=========================================================
Test script to compute  evoked data
=========================================================

"""
# Authors: Alexandre Gramfort <gramfort@nmr.mgh.harvard.edu>
#
# License: BSD (3-clause)

print __doc__

import sys
sys.path.insert(1,'/imaging/local/software/mne_python/latest_v0.11')

import numpy as np
import mne

###############################################################################
main_path = '/imaging/rf02/Semnet/'
data_path = '/imaging/rf02/Semnet/'	# where subdirs for MEG data are
#os.chdir('/home/rf02/rezvan/test1')

# get indices for subjects to be processed from command line input
# 
print sys.argv
subject_inds = []
for ss in sys.argv[1:]:
 subject_inds.append( int( ss ) )
#/imaging/rf02/Semnet/meg16_0022/160208/
#subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11]
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
# 2 4 5 6 14
bad_channels=(['EEG039','EEG043'],              #8 'EEG003',?
              ['EEG023', 'EEG034','EEG039', 'EEG041','EEG047'], #10 'EEG022',
              ['EEG003', 'EEG007', 'EEG008', 'EEG027','EEG046', 'EEG067','EEG070'],#11 keep 70?
              ['EEG020', 'EEG055'],   #12
              ['EEG044', 'EEG045','EEG055', 'EEG057', 'EEG059', 'EEG060'],#13
              ['EEG038', 'EEG039','EEG073'],    #15
              ['EEG044','EEG045'],    #16
              ['EEG002', 'EEG045','EEG046'],    #17
              ['EEG029','EEG039','EEG067'],    #18 27?
              ['EEG033','EEG034', 'EEG044','EEG045','EEG046'],    #19
              ['EEG039','EEG045'],    #20
              [''],    #21
              [''],    #22 not checked yet
              ['EEG033'], #24


                
)
ll = []
for ss in subject_inds:
 ll.append(list_all[ss])

print "ll:"
print ll

stim_delay = 0.034 # delay in s
##loop across subjects...
for ii, meg in enumerate(ll):
	
	print subject_inds[0]

## Get raw data
	raw_fname = main_path + meg + 'block_localisers_tsss_ds_raw.fif'
	raw = mne.io.Raw(raw_fname, preload=True)#, preload=True
## interpolating bad channels
#	raw.info['bads']=bad_channels[subject_inds[0]]
# 
#	if raw.info['bads']:
#         raw.interpolate_bads(reset_bads=True)
	print "raw loaded"
	
	include = []
	exclude = raw.info['bads']#[]; raw.info['bads']=[] # bads

	picks = mne.pick_types(raw.info, meg=True, eeg=True, eog=False,stim=False, include=include, exclude=exclude)
	print "picks"
	#raw.notch_filter(freqs=np.arange(50,251,50), picks=picks, method='fft', filter_length='50s',trans_bandwidth=0.1)
	raw.filter(l_freq=0.1, h_freq=48,  l_trans_bandwidth=0.05, picks=picks, method='fft',filter_length='200s')#l_trans_bandwidth=0.1,
	
# Writing events
	print "writing filtred"
	out_name=data_path + meg + 'block_localisers_tsss_filt_pnt1_48_raw.fif'
#cd out_path
	
	raw.save(out_name, overwrite=True)
	
