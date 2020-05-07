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
data_path = '/imaging/rf02/Semnet/'	# where subdirs for MEG data are
inv_path = '/imaging/rf02/Semnet/'
inv_fname = 'InvOp_ico5newreg_fft_pnt1_30_clean_ica_EMEG-inv.fif'#'InverseOperator_fft_1_48_clean_ica_EMEG-inv.fif'
event_path = '/imaging/rf02/Semnet/'
print sys.argv
subject_inds = []
for ss in sys.argv[1:]:
 subject_inds.append( int( ss ) )

subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, 17, 18]
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
            '/meg16_0125/160712/', #24 
            ]
# subjects names used for MRI data
subjects=['MRI_meg16_0030' ,#0
            'MRI_meg16_0032' ,#1
            'MRI_meg16_0034' ,#2
            'MRI_meg16_0035' ,#3
            'MRI_meg16_0042' ,#4
            'MRI_meg16_0045' ,#5
            'MRI_meg16_0052' ,#6
            'MRI_meg16_0056' ,#7
    	       'MRI_meg16_0069' ,#8
 	       'MRI_meg16_0070' ,#9
            'MRI_meg16_0072' ,#10
            'MRI_meg16_0073' ,#11
            'MRI_meg16_0075' ,#12
            'MRI_meg16_0078' ,#13
            'MRI_meg16_0082' ,#14
            'MRI_meg16_0086' ,#15
            'MRI_meg16_0097' ,#16
            'MRI_meg16_0122' ,#17
            'MRI_meg16_0125' ,#18
            ]

ll = []
for ss in subject_inds:
 ll.append(list_all[ss])

print "ll:"
print ll
n_subjects=19
stim_delay = 0.034 # delay in s ((NOTE! not in the example, taken from Olaf's script. Rezvan))
tmin, tmax = -0.3, 0.6
#evoked_sl=range(19)
#
##SemLoc
#for ii, meg in enumerate(ll):
#    print ii
#    #cov_path=data_path+meg+'noise_cov'
#    #noise_cov=mne.read_cov(cov_path)
#    raw_fname = data_path + meg + 'SemDec_blocks_tsss_filt_pnt1_30_ica_raw.fif' #clean_ssp_
#    event_fname = event_path + meg + 'SemDec_blocks_tsss_filt_pnt1_30_ica_raw-eve.fif'
#    subjects_dir = data_path 
#    
#    tmin = -0.3
#    tmax = 0.6  # Use a lower tmax to reduce multiple comparisons
#    #   Setup for reading the raw data
#    raw = io.Raw(raw_fname)
#    events = mne.read_events(event_fname)
#    stim_delay = 0.034 # delay in s
#    events[:,0] = events[:,0] +np.round( raw.info['sfreq']*stim_delay )   
#    ###############################################################################
#    # Read epochs for all channels, removing a bad one
#    print "epochs"
#    picks = mne.pick_types(raw.info, eeg=True, meg=True, eog=False, exclude='bads')
#    if ii in np.array([3,5]):
#        reject = dict(eeg=120e-6, grad=200e-12, mag=4e-12) #eog=150e-6,
#    else:
#        reject = dict(eeg=100e-6, grad=200e-12, mag=4e-12) #eog=150e-6,
#    
#    event_id={'visual': 1, 'hear': 2, 'hand': 3}#, 'pwordc': 6}
#    epochs = mne.Epochs(raw, events, event_id, tmin, tmax, picks=picks,baseline=(None, 0), proj=True, reject=reject, preload=True)
#    evoked_sl[ii]=epochs.average() 
#    print len(epochs['visual']); print len(epochs['hear']);print len(epochs['hand']);# print len(epochs['pwordc'])
#    perc=(len(epochs['visual'])+len(epochs['hear'])+ len(epochs['hand'])+ len(epochs['pwordc']))/600.
#    print perc
out_path_sl=data_path+'evoked_listall_semnet_pnt1_30_m300_600_words.npy'
#esl_ar=np.asarray(evoked_sl)
#np.save(out_path_sl,esl_ar)
esl_ar=np.load(out_path_sl)

evoked_sl=list(esl_ar)

snr_rms_slmn=np.zeros((901,19))
snr_rms_slstd=np.zeros((901,19))

for cntr in range(len(evoked_sl)):
    
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

xvals=np.arange(-300,601,1)
import matplotlib.pyplot as plt
#plt.hold(True) 
plt.plot(xvals,snr_rms_slstd[:,:].mean(axis=1),color='black')#,label='Only Concrete Cats')
plt.xlabel('time (ms)')
plt.ylabel('RMS of SNR')
plt.legend()
plt.show()
