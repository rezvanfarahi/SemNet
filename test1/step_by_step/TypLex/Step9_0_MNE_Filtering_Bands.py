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
sys.path.insert(1,'/imaging/local/software/python_packages/scikit-learn/v0.13.1/lib.linux-x86_64-2.7/sklearn/externals')
import joblib
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

subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
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

## Get raw data
	raw_fname = data_path + meg + 'typlex_pw_ssstf_fft_1_48_clean_ica_raw.fif'
	raw = mne.io.Raw(raw_fname, preload=True)#, preload=True
## interpolating bad channels
	if raw.info['bads']:
         raw.interpolate_bads_eeg()
	print "raw loaded"
	
	include = []
	exclude = []; raw.info['bads']=[] # bads
	raw_theta=raw.copy()
	raw_alpha=raw.copy()
	raw_beta=raw.copy()
	raw_gamma=raw.copy()


	picks = mne.pick_types(raw.info, meg=True, eeg=True, eog=True,
                                stim=False, include=include, exclude=exclude)
	print "picks"
	#raw.notch_filter(freqs=np.arange(50,251,50), picks=picks, method='fft', filter_length='50s',trans_bandwidth=0.1)
	raw_theta.filter(l_freq=4., h_freq=7., l_trans_bandwidth=0.1, picks=picks, method='fft',filter_length='40s',n_jobs=4)
	print "theta done!"
	raw_alpha.filter(l_freq=8., h_freq=12., l_trans_bandwidth=0.1, picks=picks, method='fft',filter_length='40s',n_jobs=4)
	print "alpha done!"
	raw_beta.filter(l_freq=13., h_freq=30., l_trans_bandwidth=0.1, picks=picks, method='fft',filter_length='40s',n_jobs=4)
	print "beta done!"
	raw_gamma.filter(l_freq=31., h_freq=45., l_trans_bandwidth=0.1, picks=picks, method='fft',filter_length='40s',n_jobs=4)
	print "gamma done!"

	
# Writing events
	print "writing filtred"
	out_namet=data_path + meg + 'typlex_pw_ssstf_fft_ica_theta_raw.fif'
	out_namea=data_path + meg + 'typlex_pw_ssstf_fft_ica_alpha_raw.fif'
	out_nameb=data_path + meg + 'typlex_pw_ssstf_fft_ica_beta_raw.fif'
	out_nameg=data_path + meg + 'typlex_pw_ssstf_fft_ica_gamma_raw.fif'

#cd out_path
	
	raw_theta.save(out_namet, overwrite=True)
	raw_alpha.save(out_namea, overwrite=True)
	raw_beta.save(out_nameb, overwrite=True)
	raw_gamma.save(out_nameg, overwrite=True)
	
