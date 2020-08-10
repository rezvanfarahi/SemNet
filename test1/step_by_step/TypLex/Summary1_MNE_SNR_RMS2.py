"""
=================================================================
Permutation t-test on source data with spatio-temporal clustering
=================================================================

Tests if the evoked response is significantly different between
conditions across subjects (simulated here using one subject's data).
The multiple comparisons problem is addressed with a cluster-level
permutation test across space and time.

"""

# Authors: Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
#          Eric Larson <larson.eric.d@gmail.com>
# License: BSD (3-clause)

print(__doc__)
# Russell's addition
import sys
sys.path.insert(1,'/imaging/local/software/mne_python/latest_v0.9full')
#sys.path.insert(1,'/imaging/rf02/TypLexMEG/mne_python_v8')


# for qsub
# (add ! here if needed) /imaging/local/software/anaconda/latest/x86_64/bin/python
sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
###


import numpy as np
import mne
from mne import io
from mne.minimum_norm import apply_inverse, read_inverse_operator
from mne.epochs import equalize_epoch_counts

###############################################################################
# Set parameters
data_path = '/imaging/rf02/TypLexMEG/'	# where subdirs for MEG data are
inv_path = '/imaging/rf02/TypLexMEG/'
inv_fname = 'typlex_InverseOperator_fft_pnt1_48_clean_ica_EMEG-inv.fif'#'InverseOperator_fft_1_48_clean_ica_EMEG-inv.fif'
event_path = '/group/erp/data/olaf.hauk/MEG/TypLexM/data/MF22_new/'
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
n_subjects=17
stim_delay = 0.034 # delay in s ((NOTE! not in the example, taken from Olaf's script. Rezvan))
tmin, tmax = -0.3, 0.55
X=np.zeros((20484,35, 17,2))
evoked_sl=range(17)
evoked_tl=range(17)

for ii, meg in enumerate(ll):
    print ii
    #cov_path=data_path+meg+'noise_cov'
    #noise_cov=mne.read_cov(cov_path)
    raw_fname = data_path + meg + 'typlex_pw_ssstf_fft_pnt3_48_clean_ica_raw.fif' #clean_ssp_
    event_fname = event_path + meg + 'typlex_pw_raw_ssstnew.txt'
    subjects_dir = data_path 
    events = mne.read_events(event_fname)
    if ii < 3:
	    	word_ids=np.array([2,1]); print events[events[:,2]==1,2].shape+events[events[:,2]==2,2].shape
	    	nonword_ids=np.array([6,5]); print events[events[:,2]==5,2].shape+events[events[:,2]==6,2].shape
	    	events=mne.merge_events(events,word_ids,601,replace_events=True)
	    	events=mne.merge_events(events,nonword_ids,602,replace_events=True)
	    	event_ids = {'Words': 601, 'PWords': 602}

    else:
	    	word_ids=np.array([8,7]); print events[events[:,2]==7,2].shape+events[events[:,2]==8,2].shape
	    	nonword_ids=np.array([6,5]); print events[events[:,2]==5,2].shape+events[events[:,2]==6,2].shape
	    	events=mne.merge_events(events,word_ids,601,replace_events=True)
	    	events=mne.merge_events(events,nonword_ids,602,replace_events=True)
	    	event_ids = {'Words': 601, 'PWords': 602}
    tmin = -0.3
    tmax = 0.55  # Use a lower tmax to reduce multiple comparisons
    #   Setup for reading the raw data
    raw = io.Raw(raw_fname)
    
    stim_delay = 0.034 # delay in s
    events[:,0] += np.round( raw.info['sfreq']*stim_delay )   
    ###############################################################################
    # Read epochs for all channels, removing a bad one
    print "epochs"
    picks = mne.pick_types(raw.info, eeg=True, meg=True, eog=False, exclude='bads')
#    event_id = 601  # word
    if ii in np.array([3,4,5,13]):
        reject = dict(eeg=150e-6, grad=200e-12, mag=4e-12) #eog=150e-6,
    else:
        reject = dict(eeg=120e-6, grad=200e-12, mag=4e-12) #eog=150e-6,
#    epochs1 = mne.Epochs(raw, events, event_id, tmin, tmax, picks=picks,
#    		             baseline=(None, 0), proj=True, reject=reject, preload=True)
#    
#    event_id = 602  # nonword
#    epochs2 = mne.Epochs(raw, events, event_id, tmin, tmax, picks=picks, baseline=(None, 0), proj=True, reject=reject, preload=True)
    epochs = mne.Epochs(raw, events, event_ids, tmin, tmax, picks=picks,baseline=(None, 0), proj=True, reject=reject, preload=True)
    print len(epochs['Words']); print len(epochs['PWords'])
#    equalize_epoch_counts(epochs,method='mintime')

                   
                   
                   
	#    Equalize trial counts to eliminate bias (which would otherwise be
	#    introduced by the abs() performed below)

	###############################################################################
    print "Transform to source space"
    
    fname_inv = inv_path + meg + inv_fname
    
    method = "MNE"  # use dSPM method (could also be MNE or sLORETA)
    #inverse_operator = read_inverse_operator(fname_inv)
    #sample_vertices = [s['vertno'] for s in inverse_operator['src']]
    #    Let's average and compute inverse, resampling to speed things up
    evoked_tl[ii]=epochs.average() 
    
    
##SemLoc
for ii, meg in enumerate(ll):
    print ii
    #cov_path=data_path+meg+'noise_cov'
    #noise_cov=mne.read_cov(cov_path)
    raw_fname = data_path + meg + 'semloc_ssstf_fft_pnt3_48_clean_ica_raw.fif' #clean_ssp_
    event_fname = event_path + meg + 'semloc_raw_ssstnew.txt'
    subjects_dir = data_path 
    
    tmin = -0.3
    tmax = 0.55  # Use a lower tmax to reduce multiple comparisons
    #   Setup for reading the raw data
    raw = io.Raw(raw_fname)
    events = mne.read_events(event_fname)
    stim_delay = 0.034 # delay in s
    events[:,0] += np.round( raw.info['sfreq']*stim_delay )   
    ###############################################################################
    # Read epochs for all channels, removing a bad one
    print "epochs"
    picks = mne.pick_types(raw.info, eeg=True, meg=True, eog=False, exclude='bads')
    if ii in np.array([2,3]):
        reject = dict(eeg=150e-6, grad=200e-12, mag=4e-12) #eog=150e-6,
    else:
        reject = dict(eeg=120e-6, grad=200e-12, mag=4e-12) #eog=150e-6,
    
    event_id={'cncrt_wrd': 1, 'abs_wrd': 2}
    epochs = mne.Epochs(raw, events, event_id, tmin, tmax, picks=picks,baseline=(None, 0), proj=True, reject=reject, preload=True)
    evoked_sl[ii]=epochs.average() 
    print len(epochs['cncrt_wrd']); print len(epochs['abs_wrd'])
out_path_sl=data_path+'evoked_listall_sl_pnt3_m300_550.npy'
esl_ar=np.asarray(evoked_sl)
np.save(out_path_sl,esl_ar)
esl_ar=np.load(out_path_sl)
out_path_tl=data_path+'evoked_listall_tl_pnt3_m300_550.npy'
etl_ar=np.asarray(evoked_tl)
np.save(out_path_tl,etl_ar)
etl_ar=np.load(out_path_tl)
evoked_tl=list(etl_ar)
evoked_sl=list(esl_ar)
snr_rms_tlmn=np.zeros((851,17))
snr_rms_tlstd=np.zeros((851,17))
snr_rms_slmn=np.zeros((851,17))
snr_rms_slstd=np.zeros((851,17))

for cntr in range(len(evoked_tl)):
    evoked_tl_cntr=evoked_tl[cntr]#mne.grand_average(evoked_tl)#,interpolate_bads=False
    epdatatl=evoked_tl_cntr.data 
    emtl=epdatatl.copy()#np.abs(epdata)
    #emax=np.mean(em[:,:,600:900],axis=2)#
    signaltl=emtl#[:,100:1100]#.max(axis=2)
    #bm=np.sqrt(np.mean(epdata[:,:,100:500]**2,axis=2))
    noisetl=np.mean(emtl[:,:300],axis=1)#np.mean(emtl[:,100:500],axis=1)
    noise_stdtl=np.std(emtl[:,:300],axis=1)#np.std(emtl[:,100:500],axis=1)
    
    edvtl=signaltl.copy()
    edv_stdtl=signaltl.copy()
    
    for ii in range(signaltl.shape[1]):
        edvtl[:,ii]=(signaltl[:,ii]/noisetl)**2
        edv_stdtl[:,ii]=(signaltl[:,ii]/noise_stdtl)**2
    edivtl=np.sqrt(edvtl.mean(axis=0))
    ediv_stdtl=np.sqrt(edv_stdtl.mean(axis=0))#edv_stdtl.mean(axis=0)#
    snr_rms_tlmn[:,cntr]=edivtl.copy()
    snr_rms_tlstd[:,cntr]=ediv_stdtl.copy()
    
    evoked_sl_cntr=evoked_sl[cntr]#evoked_sl_ga=mne.grand_average(evoked_sl)#,interpolate_bads=False
    epdatasl=evoked_sl_cntr.data 
    emsl=epdatasl.copy()#np.abs(epdata)
    #emax=np.mean(em[:,:,600:900],axis=2)#
    signalsl=emsl#[:,100:1100]#.max(axis=2)
    #bm=np.sqrt(np.mean(epdata[:,:,100:500]**2,axis=2))
    noisesl=np.mean(emsl[:,:300],axis=1)#np.mean(emsl[:,100:500],axis=1)
    noise_stdsl=np.std(emsl[:,:300],axis=1)#np.std(emsl[:,100:500],axis=1)
    
    edvsl=signalsl.copy()
    edv_stdsl=signalsl.copy()
    
    for ii in range(signalsl.shape[1]):
        edvsl[:,ii]=(signalsl[:,ii]/noisesl)**2
        edv_stdsl[:,ii]=(signalsl[:,ii]/noise_stdsl)**2
    edivsl=np.sqrt(edvsl.mean(axis=0))
    ediv_stdsl=np.sqrt(edv_stdsl.mean(axis=0))#edv_stdsl.mean(axis=0)#np.sqrt(edv_stdsl.mean(axis=0))
    snr_rms_slmn[:,cntr]=edivsl.copy()
    snr_rms_slstd[:,cntr]=ediv_stdsl.copy()


#snr = 3.0#np.mean(ediv)
#xvals=np.arange(-200,600,1)
#import matplotlib.pyplot as plt
#plt.plot(xvals,edivtl[200:],'r',xvals,edivsl[200:],'k')    
#plt.xlabel('time (ms)')
#plt.ylabel('RMS of SNR')
#plt.show()

xvals=np.arange(-300,551,1)
import matplotlib.pyplot as plt
plt.plot(xvals,snr_rms_tlstd[:,:].mean(axis=1),color='red',label='Lexical Decision')
plt.hold(True) 
plt.plot(xvals,snr_rms_slstd[:,:].mean(axis=1),color='black',label='Semantic Decision')
plt.xlabel('time (ms)')
plt.ylabel('RMS of SNR')
plt.legend()
plt.show()
