# Author: Martin Luessi <mluessi@nmr.mgh.harvard.edu>
#
# License: BSD (3-clause)

print(__doc__)
print "hi! this script is for spectral functional connectivity. Very nice maps are fruits which worth waiting for I bet!"
import sys

sys.path.append('/imaging/local/software/python_packages/nibabel/1.3.0')
sys.path.append('/imaging/local/software/python_packages/pysurfer/v0.4')
#sys.path.append('/imaging/local/software/python_packages/scikit-learn/v0.14.1')
sys.path.append('/imaging/local/software/mne_python/latest')

# for qsub
# (add ! here if needed) 
#sys.path.append('/imaging/local/software/EPD/latest/x86_64/bin/python')
#sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
###

import os
import numpy as np
import mne
from mne.connectivity import seed_target_indices
data_path = '/imaging/rf02/TypLexMEG/' # root directory for your MEG data
os.chdir(data_path)
event_path = '/group/erp/data/olaf.hauk/MEG/TypLexM/data/MF22_new/'  # where event files are
from scipy.stats import pearsonr
label_path = '/imaging/rf02/TypLexMEG/fsaverage/label/'
inv_path = '/imaging/rf02/TypLexMEG/'
subjects_dir = data_path 
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

ll = []
for ss in subject_inds:
 ll.append(list_all[ss])

print "ll:"
print ll




#label_name_lh = ['atlleft-lh', 'atlright-rh']


# Labels/ROIs
labellist = ['leftATL-lh', 'rightATL-rh']#, 'medialtempright-rh','medialtempleft-lh']

ii=-1
for cntt, meg in enumerate(ll):
	print cntt
	
	sct_in = data_path + meg + 'firstMorphed_SemLoc_icaclean_Concrete_Source_Evoked_m500_700'
	stcs = mne.read_source_estimate(sct_in)

	Fs=np.round(1/(stcs.times[1]-stcs.times[0]))
	stcs_d=stcs.copy().data; stcs_d=np.float64(stcs_d)
	fmin=4; fmax=8
	stcs_d=mne.filter.band_pass_filter(stcs_d,Fs=Fs,Fp1=fmin, Fp2=fmax, method='fft')
	stcs = mne.SourceEstimate(stcs_d, vertices=stcs.vertno, tmin=stcs.tmin, tstep=stcs.tstep, subject='fsaverage')
	band='theta';print "filtering theta finished"
# Now we are ready to compute the Coherence in the alpha and beta band.
# fmin and fmax specify the lower and upper freq. for each band, resp.
		
	for label_name in labellist:
		print label_name
		# frequency bands with variable number of cycles for wavelets

		fname_label = label_path  + label_name + '.label'

		label1 = mne.read_label(fname_label)
		
# Restrict the source estimate to the label in the left auditory cortex
		stc=stcs.copy()
		stc_label = stc.in_label(label1)

# Find number and index of vertex with most power
		src_pow = np.sum(stc_label.data ** 2, axis=1)
		if label_name[-2] is 'l':

			seed_vertno = stc_label.vertno[0][np.argmax(src_pow)]
			seed_idx = np.searchsorted(stc.vertno[0], seed_vertno)  # index in original stc
			print seed_idx

		else:
			print label_name
			seed_vertno = stc_label.vertno[1][np.argmax(src_pow)]
			seed_idx = stc.vertno[0].shape[0]+ np.searchsorted(stc.vertno[1], seed_vertno)  # index in original stc
			print seed_idx

# Generate index parameter for seed-based connectivity analysis
		n_sources = stc.data.shape[0]
		indices = seed_target_indices([seed_idx], np.arange(n_sources))
		tar_ind=indices[0][0]

		print "temporal connectivity correlation"
		b2=0.25
		Rc=np.zeros((stcs.data.shape[0],2))
		pc=np.zeros((stcs.data.shape[0],2))
		
		for cc in range(2):
			b1=b2-0.1+0.001; b2=b1+0.2-0.001
			datamatx=stcs.copy().crop(b1,b2).data	
			#dh=scipy.signal.hilbert(datamatx)
			#datamatx=np.abs(dh)
			for dcnt in range(datamatx.shape[0]):
				Rc[dcnt,cc],pc[dcnt,cc]=pearsonr(datamatx[dcnt,:],datamatx[tar_ind,:])
		tmin1 = 1#np.mean(freqs[0])
		tstep1 = 1#np.mean(freqs[1]) - tmin
		corr_stc = mne.SourceEstimate(Rc, vertices=stcs.vertno, tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
		data_out = data_path + meg +  'firstMorphed_SemLoc_icaclean_Concrete_Correlation_' + band + '_150_450_200ms_' + label_name[0:-3]
		
		corr_stc.save(data_out)
		p_stc = mne.SourceEstimate(Rc, vertices=stcs.vertno, tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
		data_out = data_path + meg +  'Morphed_SemLoc_icaclean_Concrete_Correlation_' + band + '_pvalues_150_450_200ms_' + label_name[0:-3]
		p_stc.save(data_out)
############
print "theta finished"	


ii=-1
for cntt, meg in enumerate(ll):
	print cntt
	
	sct_in = data_path + meg + 'firstMorphed_SemLoc_icaclean_Concrete_Source_Evoked_m500_700'
	stcs = mne.read_source_estimate(sct_in)

	Fs=np.round(1/(stcs.times[1]-stcs.times[0]))
	stcs_d=stcs.copy().data; stcs_d=np.float64(stcs_d)
	fmin=8; fmax=12
	stcs_d=mne.filter.band_pass_filter(stcs_d,Fs=Fs,Fp1=fmin, Fp2=fmax, l_trans_bandwidth=0.2, h_trans_bandwidth=0.2, method='fft')
	stcs = mne.SourceEstimate(stcs_d, vertices=stcs.vertno, tmin=stcs.tmin, tstep=stcs.tstep, subject='fsaverage')
	band='alpha';print "filtering alpha finished"
# Now we are ready to compute the Coherence in the alpha and beta band.
# fmin and fmax specify the lower and upper freq. for each band, resp.
		
	for label_name in labellist:
		print label_name
		# frequency bands with variable number of cycles for wavelets

		fname_label = label_path  + label_name + '.label'

		label1 = mne.read_label(fname_label)
		
# Restrict the source estimate to the label in the left auditory cortex
		stc=stcs.copy()
		stc_label = stc.in_label(label1)

# Find number and index of vertex with most power
		src_pow = np.sum(stc_label.data ** 2, axis=1)
		if label_name[-2] is 'l':

			seed_vertno = stc_label.vertno[0][np.argmax(src_pow)]
			seed_idx = np.searchsorted(stc.vertno[0], seed_vertno)  # index in original stc
			print seed_idx

		else:
			print label_name
			seed_vertno = stc_label.vertno[1][np.argmax(src_pow)]
			seed_idx = stc.vertno[0].shape[0]+ np.searchsorted(stc.vertno[1], seed_vertno)  # index in original stc
			print seed_idx

# Generate index parameter for seed-based connectivity analysis
		n_sources = stc.data.shape[0]
		indices = seed_target_indices([seed_idx], np.arange(n_sources))
		tar_ind=indices[0][0]

		print "temporal connectivity correlation"
		b2=0.25
		Rc=np.zeros((stcs.data.shape[0],2))
		pc=np.zeros((stcs.data.shape[0],2))
		
		for cc in range(2):
			b1=b2-0.1+0.001; b2=b1+0.2-0.001
			datamatx=stcs.copy().crop(b1,b2).data	
			#dh=scipy.signal.hilbert(datamatx)
			#datamatx=np.abs(dh)
			for dcnt in range(datamatx.shape[0]):
				Rc[dcnt,cc],pc[dcnt,cc]=pearsonr(datamatx[dcnt,:],datamatx[tar_ind,:])
		tmin1 = 1#np.mean(freqs[0])
		tstep1 = 1#np.mean(freqs[1]) - tmin
		corr_stc = mne.SourceEstimate(Rc, vertices=stcs.vertno, tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
		data_out = data_path + meg +  'firstMorphed_SemLoc_icaclean_Concrete_Correlation_' + band + '_150_450_200ms_' + label_name[0:-3]
		
		corr_stc.save(data_out)
		p_stc = mne.SourceEstimate(Rc, vertices=stcs.vertno, tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
		data_out = data_path + meg +  'Morphed_SemLoc_icaclean_Concrete_Correlation_' + band + '_pvalues_150_450_200ms_' + label_name[0:-3]
		p_stc.save(data_out)
############
print "alpha finished"	


ii=-1
for cntt, meg in enumerate(ll):
	print cntt
	
	sct_in = data_path + meg + 'firstMorphed_SemLoc_icaclean_Concrete_Source_Evoked_m500_700'
	stcs = mne.read_source_estimate(sct_in)

	Fs=np.round(1/(stcs.times[1]-stcs.times[0]))
	stcs_d=stcs.copy().data; stcs_d=np.float64(stcs_d)
	fmin=13; fmax=30
	stcs_d=mne.filter.band_pass_filter(stcs_d,Fs=Fs,Fp1=fmin, Fp2=fmax, l_trans_bandwidth=0.2, h_trans_bandwidth=0.2, method='fft')
	stcs = mne.SourceEstimate(stcs_d, vertices=stcs.vertno, tmin=stcs.tmin, tstep=stcs.tstep, subject='fsaverage')
	band='beta';print "filtering beta finished"
# Now we are ready to compute the Coherence in the alpha and beta band.
# fmin and fmax specify the lower and upper freq. for each band, resp.
		
	for label_name in labellist:
		print label_name
		# frequency bands with variable number of cycles for wavelets

		fname_label = label_path  + label_name + '.label'

		label1 = mne.read_label(fname_label)
		
# Restrict the source estimate to the label in the left auditory cortex
		stc=stcs.copy()
		stc_label = stc.in_label(label1)

# Find number and index of vertex with most power
		src_pow = np.sum(stc_label.data ** 2, axis=1)
		if label_name[-2] is 'l':

			seed_vertno = stc_label.vertno[0][np.argmax(src_pow)]
			seed_idx = np.searchsorted(stc.vertno[0], seed_vertno)  # index in original stc
			print seed_idx

		else:
			print label_name
			seed_vertno = stc_label.vertno[1][np.argmax(src_pow)]
			seed_idx = stc.vertno[0].shape[0]+ np.searchsorted(stc.vertno[1], seed_vertno)  # index in original stc
			print seed_idx

# Generate index parameter for seed-based connectivity analysis
		n_sources = stc.data.shape[0]
		indices = seed_target_indices([seed_idx], np.arange(n_sources))
		tar_ind=indices[0][0]

		print "temporal connectivity correlation"
		b2=0.25
		Rc=np.zeros((stcs.data.shape[0],2))
		pc=np.zeros((stcs.data.shape[0],2))
		
		for cc in range(2):
			b1=b2-0.1+0.001; b2=b1+0.2-0.001
			datamatx=stcs.copy().crop(b1,b2).data	
			#dh=scipy.signal.hilbert(datamatx)
			#datamatx=np.abs(dh)
			for dcnt in range(datamatx.shape[0]):
				Rc[dcnt,cc],pc[dcnt,cc]=pearsonr(datamatx[dcnt,:],datamatx[tar_ind,:])
		tmin1 = 1#np.mean(freqs[0])
		tstep1 = 1#np.mean(freqs[1]) - tmin
		corr_stc = mne.SourceEstimate(Rc, vertices=stcs.vertno, tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
		data_out = data_path + meg +  'firstMorphed_SemLoc_icaclean_Concrete_Correlation_' + band + '_150_450_200ms_' + label_name[0:-3]
		
		corr_stc.save(data_out)
		p_stc = mne.SourceEstimate(Rc, vertices=stcs.vertno, tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
		data_out = data_path + meg +  'Morphed_SemLoc_icaclean_Concrete_Correlation_' + band + '_pvalues_150_450_200ms_' + label_name[0:-3]
		p_stc.save(data_out)
############
print "beta finished"	

ii=-1
for cntt, meg in enumerate(ll):
	print cntt
	
	sct_in = data_path + meg + 'firstMorphed_SemLoc_icaclean_Concrete_Source_Evoked_m500_700'
	stcs = mne.read_source_estimate(sct_in)

	Fs=np.round(1/(stcs.times[1]-stcs.times[0]))
	stcs_d=stcs.copy().data; stcs_d=np.float64(stcs_d)
	fmin=31; fmax=45
	stcs_d=mne.filter.band_pass_filter(stcs_d,Fs=Fs,Fp1=fmin, Fp2=fmax,l_trans_bandwidth=0.2, h_trans_bandwidth=0.2, method='fft')
	stcs = mne.SourceEstimate(stcs_d, vertices=stcs.vertno, tmin=stcs.tmin, tstep=stcs.tstep, subject='fsaverage')
	band='gamma';print "filtering gamma finished"
# Now we are ready to compute the Coherence in the alpha and beta band.
# fmin and fmax specify the lower and upper freq. for each band, resp.
		
	for label_name in labellist:
		print label_name
		# frequency bands with variable number of cycles for wavelets

		fname_label = label_path  + label_name + '.label'

		label1 = mne.read_label(fname_label)
		
# Restrict the source estimate to the label in the left auditory cortex
		stc=stcs.copy()
		stc_label = stc.in_label(label1)

# Find number and index of vertex with most power
		src_pow = np.sum(stc_label.data ** 2, axis=1)
		if label_name[-2] is 'l':

			seed_vertno = stc_label.vertno[0][np.argmax(src_pow)]
			seed_idx = np.searchsorted(stc.vertno[0], seed_vertno)  # index in original stc
			print seed_idx

		else:
			print label_name
			seed_vertno = stc_label.vertno[1][np.argmax(src_pow)]
			seed_idx = stc.vertno[0].shape[0]+ np.searchsorted(stc.vertno[1], seed_vertno)  # index in original stc
			print seed_idx

# Generate index parameter for seed-based connectivity analysis
		n_sources = stc.data.shape[0]
		indices = seed_target_indices([seed_idx], np.arange(n_sources))
		tar_ind=indices[0][0]

		print "temporal connectivity correlation"
		b2=0.25
		Rc=np.zeros((stcs.data.shape[0],2))
		pc=np.zeros((stcs.data.shape[0],2))
		
		for cc in range(2):
			b1=b2-0.1+0.001; b2=b1+0.2-0.001
			datamatx=stcs.copy().crop(b1,b2).data	
			#dh=scipy.signal.hilbert(datamatx)
			#datamatx=np.abs(dh)
			for dcnt in range(datamatx.shape[0]):
				Rc[dcnt,cc],pc[dcnt,cc]=pearsonr(datamatx[dcnt,:],datamatx[tar_ind,:])
		tmin1 = 1#np.mean(freqs[0])
		tstep1 = 1#np.mean(freqs[1]) - tmin
		corr_stc = mne.SourceEstimate(Rc, vertices=stcs.vertno, tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
		data_out = data_path + meg +  'firstMorphed_SemLoc_icaclean_Concrete_Correlation_' + band + '_150_450_200ms_' + label_name[0:-3]
		
		corr_stc.save(data_out)
		p_stc = mne.SourceEstimate(Rc, vertices=stcs.vertno, tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
		data_out = data_path + meg +  'Morphed_SemLoc_icaclean_Concrete_Correlation_' + band + '_pvalues_150_450_200ms_' + label_name[0:-3]
		p_stc.save(data_out)
############
print "gamma finished"	
