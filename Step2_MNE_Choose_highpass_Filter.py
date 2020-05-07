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
import scipy.io as scio
#os.chdir('/home/rf02/rezvan/test1')

# get indices for subjects to be processed from command line input
# 
print sys.argv
subject_inds = []
for ss in sys.argv[1:]:
 subject_inds.append( int( ss ) )
#/imaging/rf02/Semnet/meg16_0022/160208/
#subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]
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
            '/meg16_0122/160707/', #22 LD
            '/meg16_0125/160712/', #24 LD
            ]

bad_channels=(['EEG008','EEG028'],#0 
              ['EEG067'],#1
              ['EEG027', 'EEG028'],#3 remove?
              ['EEG013', 'EEG038', 'EEG039','EEG073'],#5
              ['EEG003', 'EEG004','EEG022', 'EEG023', 'EEG037', 'EEG038', 'EEG045', 'EEG046','EEG059', 'EEG072'],#6
              ['EEG002', 'EEG034', 'EEG045','EEG046'],    #9 17, 18, 55, 57?
              ['EEG023', 'EEG034','EEG039', 'EEG041','EEG047'], #10 'EEG022',
              ['EEG003', 'EEG007', 'EEG008', 'EEG027','EEG046', 'EEG067','EEG070'],#11 keep 70?
              ['EEG020', 'EEG055'],   #12
              ['EEG044', 'EEG045','EEG055', 'EEG057', 'EEG059', 'EEG060'],#13
              ['EEG038', 'EEG039','EEG073'],    #15
              ['EEG044','EEG045'],    #16
              ['EEG002', 'EEG045','EEG046'],    #17
              ['EEG029','EEG039','EEG067'],    #18
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
lfreq=1
stim_delay = 0.034 # delay in s
##loop across subjects...
snr_vect=np.zeros((len(list_all),))
for ii, meg in enumerate(ll):
    print subject_inds[0]

## Get raw data
    raw_fname1 = main_path + meg + 'block_fruit_tsss_raw.fif'
    raw1 = mne.io.Raw(raw_fname1, preload=True)#, preload=True
    raw_fname2 = main_path + meg + 'block_milk_tsss_raw.fif'
    raw2 = mne.io.Raw(raw_fname2, preload=True)#, preload=True
    raw_fname3 = main_path + meg + 'block_odour_tsss_raw.fif'
    raw3 = mne.io.Raw(raw_fname3, preload=True)#, preload=True
    
    print "Reading events"   
    events_fname1 = data_path + meg + 'block_fruit_tsss_raw-eve.fif'#orig_path + meg + 'semloc_raw_ssstnew.txt'
    events1 = mne.read_events(events_fname1)
    events_fname2 = data_path + meg + 'block_milk_tsss_raw-eve.fif'#orig_path + meg + 'semloc_raw_ssstnew.txt'
    events2 = mne.read_events(events_fname2)
    events_fname3 = data_path + meg + 'block_odour_tsss_raw-eve.fif'#orig_path + meg + 'semloc_raw_ssstnew.txt'
    events3 = mne.read_events(events_fname3)
    
    raws=[raw1,raw2,raw3]
    eventss=[events1,events2,events3]
    raw,events=mne.concatenate_raws(raws,events_list=eventss)
## interpolating bad channels
#	raw.info['bads']=bad_channels[subject_inds[0]]
# 
#	if raw.info['bads']:
#         raw.interpolate_bads(reset_bads=True)
    print "raw loaded"
    tmin=-0.3
    tmax=0.6
    include = []
    exclude = raw.info['bads']#[]; raw.info['bads']=[] # bads
    picks = mne.pick_types(raw.info, meg=True, eeg=True, eog=False,stim=False, include=include, exclude=exclude)
    print "picks"
	#raw.notch_filter(freqs=np.arange(50,251,50), picks=picks, method='fft', filter_length='50s',trans_bandwidth=0.1)
    raw.filter(l_freq=lfreq, h_freq=48,  l_trans_bandwidth=0.05, picks=picks, method='fft',filter_length='200s')#l_trans_bandwidth=0.1,
    picks = mne.pick_types(raw.info, eeg=True, meg=True, eog=True, exclude='bads')
    print "Epoching"
    if subject_inds[0] in np.array([3,5]):#np.array([4,8,9]):
        reject = dict(eeg=150e-6, grad=200e-12, mag=4e-12, eog=150e-6)
    else:
        reject = dict(eeg=120e-6, grad=200e-12, mag=4e-12, eog=150e-6)
    event_ids = {'visual': 1, 'hear': 2, 'hand': 3, 'neutral': 4, 'emotional': 5,'pwordc': 6}
    epochs = mne.Epochs(raw, events, event_ids, tmin, tmax, picks=picks, proj=True, baseline=(None, 0), reject=reject)
    evoked=epochs.average()    
    print evoked.nave
    epdata=evoked.data 
    em=np.abs(epdata)
    #emax=np.mean(em[:,:,600:900],axis=2)#
    sa=em[:,350:750]#.max(axis=2)
    #bm=np.sqrt(np.mean(epdata[:,:,100:500]**2,axis=2))
    bv=np.var(em[:,50:300],axis=1)
    edv=sa.copy()
    for jj in range(sa.shape[1]):
        edv[:,jj]=(sa[:,jj]**2)/bv
    ediv=np.sqrt(edv)
    snr = np.mean(ediv)
    path_mtx=data_path+meg+'filt1_snr.mat'
    scio.savemat(path_mtx,{'snr':snr})
#ltc=scipy.io.loadmat(path_mtx,mat_dtype=True); stc_mean_data=np.squeeze(ltc['stc_mean_data'])

#    snr_vect[ii]=snr
    print snr
# Writing events
    print "writing filtred"
#	out_name=data_path + meg + 'block_LD_tsss_filt_pnt1_48_raw.fif'
##cd out_path
#	
#	raw.save(out_name, overwrite=True)
#	
