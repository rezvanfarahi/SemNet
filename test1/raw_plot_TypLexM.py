print(__doc__)
import sys
sys.path.append('/imaging/local/software/python_packages/nibabel/1.3.0')
sys.path.append('/imaging/local/software/python_packages/pysurfer/v0.4')
sys.path.append('/imaging/local/software/python_packages/scikit-learn/v0.14.1')
sys.path.append('/imaging/local/software/mne_python/latest')

# for qsub
# (add ! here if needed) 
#sys.path.append('/imaging/local/software/anaconda/latest/x86_64/bin/python')
sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
###
import mne
import pylab as pl
from mne import io
from mne.io import Raw
from mne.preprocessing import read_ica
data_path = '/imaging/rf02/TypLexMEG/'
print sys.argv
subject_inds = [6]
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



# Labels/ROIs
#labellist = ['atlleft-lh']#, 'atlright-rh'
tmin, tmax = -0.5, 0.7
reject = dict(eeg=120e-6, grad=200e-12, mag=4e-12) #eog=150e-6,
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
for meg in ll:
	print meg
	fname = data_path + meg + 'semloc_ssstf_fft_1_48_raw.fif'
	fname1 = data_path + meg + 'semloc_ssstf_fft_1_48_clean_ica_raw.fif'

	raw = mne.io.Raw(fname)
	raw1 = mne.io.Raw(fname1)

	#ica_path=data_path+meg+'eeg-ica.fif'
	#ica_eeg=read_ica(ica_path)
	#ica_path=data_path+meg+'meg-ica.fif'
	#ica_meg=read_ica(ica_path)
	# Set up pick list: MEG + STI101 - bad channels
	#include = ['STI101']
	#raw.info['bads'] += ['MEG 2443', 'EEG 053']  # bad channels + 2 more
	print "loaded"
	picks = mne.pick_types(raw.info, meg=False, eeg=True, stim=False, eog=False, exclude='bads')
	#raw=ica_eeg.apply(raw)
	#print "ica eeg applied"
	#raw=ica_meg.apply(raw)
	#print "ica meg applied"

	
	#raw.filter(l_freq=1.5, h_freq=48, l_trans_bandwidth=.2, picks=picks, method='fft')
	some_picks = picks[:]  # take 5 first
	data=raw.copy()
	data, times = data[some_picks, 1000:]
	data=data.transpose();

	# save 150s of MEG data in FIF file
	#raw.save('sample_audvis_meg_raw.fif', tmin=0, tmax=150, picks=picks,
	 #       overwrite=True)

	###############################################################################
	# Show MEG data
	pl.plot(data)
	pl.show()	
 
	data1=raw1.copy()      
	data1, times1 = data1[some_picks, 1000:]
	data1=data1.transpose();

	# save 150s of MEG data in FIF file
	#raw.save('sample_audvis_meg_raw.fif', tmin=0, tmax=150, picks=picks,
	 #       overwrite=True)

	###############################################################################
	# Show MEG data
	pl.plot(data1)
	pl.show()	
	"""
	data=raw1.copy().apply_proj()
	print "ssp applied"
	data, times = data[some_picks, 2000:]
	data=data.transpose();

	# save 150s of MEG data in FIF file
	#raw.save('sample_audvis_meg_raw.fif', tmin=0, tmax=150, picks=picks,
	 #       overwrite=True)

	###############################################################################
	# Show MEG data
	pl.plot(data)
	pl.show()
	
	"""
