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
sys.path.insert(1,'/imaging/local/software/mne_python/latest_v0.9full')

import numpy as np
import mne

###############################################################################
main_path = '/group/erp/data/olaf.hauk/MEG/TypLexM/data/MF22_new/'
data_path = '/imaging/rf02/TypLexMEG/'	# where subdirs for MEG data are
#os.chdir('/home/rf02/rezvan/test1')

# get indices for subjects to be processed from command line input
# 
print sys.argv
subject_inds = []
for ss in sys.argv[1:]:
 subject_inds.append( int( ss ) )

#subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
print "subject_inds:"
print subject_inds
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
'meg11_0147/110603/'
]


ll = []
for ss in subject_inds:
 ll.append(list_all[ss])

print "ll:"
print ll

stim_delay = 0.034 # delay in s
##loop across subjects...
for ii, meg in enumerate(ll):
	
	print meg

# Get raw data
	raw_fname = main_path + meg + 'semloc_raw_ssst.fif'
	raw = mne.io.Raw(raw_fname, preload=True); print raw.info['bads']#, preload=True; 
## interpolating bad channels
	if raw.info['bads']:
         raw.interpolate_bads(reset_bads=True)
	print "raw loaded"
	
	include = []
	exclude = []; raw.info['bads']=[] # bads

	picks = mne.pick_types(raw.info, meg=True, eeg=True, eog=False,
                                stim=False, include=include, exclude=exclude)
	print "picks"
	#raw.notch_filter(freqs=np.arange(50,251,50), picks=picks, method='fft', filter_length='50s',trans_bandwidth=0.1)
	#raw.filter(l_freq=0.5, h_freq=48,  l_trans_bandwidth=0.1, picks=picks, method='fft',filter_length='50s')#raw.filter(l_freq=1, h_freq=48, l_trans_bandwidth=0.1, picks=picks, method='fft',filter_length='40s')
	#iir_params=dict(order=4, ftype='butter')
	#iir_params=mne.filter.construct_iir_filter(iir_params,f_pass=[0.3,48.0],f_stop=None,sfreq=1000.0,btype='bandpass',return_copy=False)
	#raw.filter( l_freq=0.3,h_freq=48, picks=picks, method='iir',iir_params=iir_params)#raw.filter(l_freq=1, h_freq=48, l_trans_bandwidth=0.1, picks=picks, method='fft',filter_length='40s')
 
	raw.filter( l_freq=0.3,h_freq=48.0, l_trans_bandwidth=0.05, picks=picks, method='fft',filter_length='200s')#raw.filter(l_freq=1, h_freq=48, l_trans_bandwidth=0.1, picks=picks, method='fft',filter_length='40s')
	
# 1Writing events
#	raw_fname = data_path + meg + 'semloc_ssstf_fft_1_48_clean_ica_raw.fif'
#	raw = mne.io.Raw(raw_fname, preload=True)
#	raw.resample(200)  
#	out_name=data_path + meg + 'semloc_ssstf_fft_1_48_resample200_clean_ica_raw.fif'
   
	print "writing filtred"
	out_name=data_path + meg + 'semloc_ssstf_fft_pnt3_48_raw.fif'
#cd out_path
	
	raw.save(out_name, overwrite=True)
	
