"""
Preprocessing MEG and EEG data using filters and ICA
====================================================

This examples filters MEG and EEG data and subsequently computes separate
ICA solutions for MEG and EEG. An html report is created, which includes
diagnostic plots of the preprocessing steps.
"""
# Author: Denis A. Engemann <denis.engemann@gmail.com>
# License: BSD (3-clause)
import sys
sys.path.insert(1,'/imaging/local/software/mne_python/v0.11')
#sys.path.insert(1,'/imaging/local/software/mne_python/latest_v0.15')

sys.path.insert(1,'/home/rf02/semnet/meeg-preprocessing')
sys.path.append('/home/rf02/rezvan/test1')
#sys.path.append('/imaging/local/software/python_packages/scikit-learn/v0.14.1')
#sys.path.insert(1,'/imaging/local/software/python_packages/scikit-learn/v0.16.1')
sys.path.insert(1,'/imaging/rf02/scikit-learn-0.16.1')
sys.path.insert(1,'/home/rf02/.local/lib/python2.7/site-packages')
# for qsub
# (add ! here if needed) 
#sys.path.append('/imaging/local/software/anaconda/latest/x86_64/bin/python')
#sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
###
import os.path as op
import mne
#reload(mne)
import sklearn
#import pylab as pl
#from mne import (fiff, read_evokeds, equalize_channels,read_proj, read_selection)

#from mne import filter
#import meeg_preprocessing
from preprocessing3 import compute_ica#from meeg_preprocessing.preprocessing import compute_ica

from utils import get_data_picks, setup_provenance
from mne.preprocessing import read_ica


data_path = '/imaging/rf02/Semnet/'
print sys.argv
subject_inds = []
for ss in sys.argv[1:]:
   subject_inds.append( int( ss ) )

subject_inds=[20]#0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]#0, 
print "subject_inds:"
print subject_inds
print "No rejection"

list_all =  ['/meg16_0030/160216/', #0
            '/meg16_0032/160218/', #1
            '/meg16_0033/160218/', #2
            '/meg16_0034/160219/', #3
            '/meg16_0035/160222/', #4
            '/meg16_0039/160225/', #5
            '/meg16_0041/160226/', #6
            '/meg16_0042/160229/', #7
            '/meg16_0045/160303/', #8
            '/meg16_0047/160304/', #9
            '/meg16_0052/160310/', #10
            '/meg16_0056/160314/',#11
            '/meg16_0069/160405/',#12 
            '/meg16_0070/160407/', #13
            '/meg16_0071/160407/', #14
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

subjects=['S000',
'S001',
'S002',
'S003',
'S004',
'S005',
'S006',
'S007',
'S008',
'S009',
'S010',
'S011',
'S012',
'S013',
'S014',
'S015',
'S016',
'S017',
'S018',
'S019',
'S020',
'S021',
'S022',
'S023',
'S024',
]
bad_channels=(['EEG008','EEG028'],#0 
              ['EEG067'],#1
              ['EEG004','EEG004', 'EEG027','EEG037','EEG045','EEG050', 'EEG072', 'EEG073', 'EEG074'],#2 remove?
              ['EEG027', 'EEG028'],#3 remove?
              ['EEG003', 'EEG007', 'EEG008', 'EEG027', 'EEG045', 'EEG057','EEG070'],#4
              ['EEG013', 'EEG038', 'EEG039','EEG073'],#5
              ['EEG003', 'EEG004','EEG022', 'EEG023', 'EEG037', 'EEG038', 'EEG045', 'EEG046','EEG059', 'EEG072'],#6
              ['EEG007', 'EEG008', 'EEG027', 'EEG070'],#7 also 55 56?
              ['EEG039','EEG043'],              #8 'EEG003',?
              ['EEG002', 'EEG034', 'EEG045','EEG046'],    #9 17, 18, 55, 57?
              ['EEG023', 'EEG034','EEG039', 'EEG041','EEG047'], #10 'EEG022',
              ['EEG003', 'EEG007', 'EEG008', 'EEG027','EEG046', 'EEG067','EEG070'],#11 keep 70?
              ['EEG020', 'EEG055'],   #12
              ['EEG044', 'EEG045','EEG055', 'EEG057', 'EEG059', 'EEG060'],#13
              ['EEG019', 'EEG039'],    #14
              ['EEG038', 'EEG039','EEG073'],    #15
              ['EEG044','EEG045'],    #16
              ['EEG002', 'EEG045','EEG046'],    #17
              ['EEG029','EEG039','EEG067'],    #18
              ['EEG033','EEG034', 'EEG044','EEG045','EEG046'],    #19
              ['EEG039','EEG045'],    #20
              [''],    #21
              [''],    #22 not checked yet
              ['EEG004','EEG016', 'EEG033','EEG047','EEG074'],    #23
              ['EEG033'], #24

ll = []
for ss in subject_inds:
 ll.append(list_all[ss])

print "ll:"
print ll



# Labels/ROIs
#labellist = ['atlleft-lh']#, 'atlright-rh'
tmin, tmax = -0.5, 0.7
event_ids = {'cncrt_wrd': 1, 'abs_wrd': 2}
# frequency bands with variable number of cycles for wavelets
stim_delay = 0.034 # delay in s
"""
frequencies = np.arange(8, 45, 1)
n_cycles = frequencies / float(7)
n_cycles[frequencies<=20] = 2
n_cycles[frequencies<=20] += np.arange(0,1,0.077)
    # n_cycles[np.logical_and((frequencies>10), (frequencies<=20))] = 3
"""

##loop across subjects...
for cnt, meg in enumerate(ll):
	subject=subjects[subject_inds[cnt]]
	print meg
	fname = data_path + meg + 'block_milk_tsss_filtnew_pnt1_30_raw.fif'

	raw = mne.io.Raw(fname, preload=True)
	raw_orig = raw.copy()#mne.io.Raw(fname, preload=True)
	print "raw loaded"
	include = [] ; exclude=raw.info['bads'] # or stim channels ['STI 014']
	picks = mne.pick_types(raw.info, meg=True, eeg=True, eog=False,stim=False, include=include, exclude=exclude)
	raw.filter(l_freq=1, h_freq=30,  l_trans_bandwidth=0.5,h_trans_bandwidth=2, picks=picks, method='fft',filter_length='200s')#l_trans_bandwidth=0.1,
#	raw.filter(l_freq=0.1, h_freq=30,  l_trans_bandwidth='auto',h_trans_bandwidth='auto', picks=picks, method='fft',filter_length='auto')#,phase='zero-double')
    

	################################################################################
	# jobs and runtime performance
	n_jobs = 1
	"""
	################################################################################
	# Filters

	# can be a list for repeated filtering
	filter_params = [dict(l_freq=1.0, h_freq=100, n_jobs=n_jobs, method='fft',
		              l_trans_bandwidth=0.1, h_trans_bandwidth=0.5)]

	notch_filter_params = dict(freqs=(50, 100, 150,))

	# margins for the PSD plot
	plot_fmin = 0.0
	plot_fmax = 120.0
	"""
	################################################################################
	# ICA

	n_components = 0.999
	# comment out to select ICA components via rank (useful with SSSed data):
	# n_components = 'rank'

	ica_meg_combined = True  # esimtate combined MAG and GRADs
	decim = 5  # decimation
	n_max_ecg, n_max_eog = 3, 3  # limit components detected due to ECG / EOG
	ica_reject = {'mag': 5e-12, 'grad': 5000e-13, 'eeg': 800e-6}

	################################################################################
	# Report

	img_scale = 1.5  # make big PNG images

	"""
	fig, report = check_apply_filter(raw, subject=subject,
		                         filter_params=filter_params,
		                         notch_filter_params=notch_filter_params,
		                         plot_fmin=plot_fmin, plot_fmax=plot_fmax,
		                         n_jobs=n_jobs, img_scale=img_scale)
	"""
	report, run_id, results_dir, logger = setup_provenance(script=__file__, results_dir='results')
	artifact_stats = dict()
     
	# get picks and iterate over channels
	for picks, ch_type in get_data_picks(raw, meg_combined=ica_meg_combined):
         if 1>0:#ch_type is 'eeg': #
             ica, _ = compute_ica(raw, picks=picks,subject=subject, n_components=n_components,n_max_ecg=n_max_ecg, n_max_eog=n_max_eog,ecg_tmin=-0.5, ecg_tmax=0.5, eog_tmin=-0.5,eog_tmax=0.5, reject=ica_reject,random_state=42,artifact_stats=artifact_stats, decim=decim, report=report, show=False, img_scale=img_scale)
             out_path=data_path+meg+'newfilt_milk_{}-ica.fif'.format(ch_type)
             print ch_type + " finished"
             ica.save(out_path)
             #raw.append(ica)
             #out_path2=data_path+meg+'semloc_ssstf_fft_1_48_clean_ica_raw.fif'
             #raw.save(out_path2)
             
             #report.save(data_path+meg+'preprocessing-report-{}.html'.format(subject), open_browser=True, overwrite=True)
	for k, v in artifact_stats.items():
         print(k, v)    
	report.save(op.join(results_dir, run_id,'preprocessing_newfilt_milk-report-{}.html'.format(subject)),open_browser=True, overwrite=True)
     
	ica_path=data_path+meg+'newfilt_milk_meg-ica.fif'
	ica=read_ica(ica_path)
	raw_orig=ica.apply(raw_orig)
	ica_path=data_path+meg+'newfilt_milk_eeg-ica.fif'
	ica=read_ica(ica_path)
	raw_orig=ica.apply(raw_orig)
	data_out_path=data_path+meg+'block_milk_tsss_filtnew_pnt1_30_ica_raw.fif'
	raw_orig.save(data_out_path, overwrite=True)


#picks_eog = mne.pick_types(raw.info, meg=True, eeg=True, ecg=False,eog=False)
#eog_channels = mne.pick_types(raw.info, meg=False, eeg=False, eog=True,ecg=False, stim=False)
#eog_epochs = mne.preprocessing.create_eog_epochs(raw, tmin=-0.2, tmax=0.5)#, picks=picks_eog, reject=ica_reject,  l_freq=1, h_freq=10)
