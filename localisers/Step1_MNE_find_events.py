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
main_path = '/imaging/rf02/Semnet/'
data_path = '/imaging/rf02/Semnet/'	# where subdirs for MEG data are
#os.chdir('/home/rf02/rezvan/test1')

# get indices for subjects to be processed from command line input
# 
print sys.argv
subject_inds = []#[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21]
for ss in sys.argv[1:]:
 subject_inds.append( int( ss ) )
#/imaging/rf02/Semnet/meg16_0022/160208/
subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13]#,1,2,3,4,5,6,7,8,9,
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
#loop over subjects...
for ii, meg in enumerate(ll):
	
     print meg

     ## Ge  t raw data
     raw_fname = main_path + meg + 'block_localisers_raw_ds.fif'#'block_odour_raw.fif'
     raw = mne.io.Raw(raw_fname, preload=True)#, preload=True
     ## interpolating bad channels
     raw.info['bads']=bad_channels[subject_inds[ii]]
     if raw.info['bads']:
         raw.interpolate_bads(reset_bads=True)
     print "raw loaded"
	
     include = []
     exclude = []; raw.info['bads']=[] # bads
     picks = mne.pick_types(raw.info, meg=True, eeg=True, eog=False,
                                stim=False, include=include, exclude=exclude)
     print "picks"
     raw_fname_new = main_path + meg + 'block_localisers_tsss_ds_raw.fif'
     events = mne.find_events(raw, stim_channel='STI101',min_duration=0.001,shortest_event=1)    
#     for evcnt in range(events.shape[0]-2):
#         if events[evcnt,2] in np.array([1,2,3,4,5,6]):
#             if events[evcnt+2,2]>1000:
#                 events[evcnt,2]=230
     mne.write_events(raw_fname_new[:-4]+'-eve.fif', events)
	
#     for ind, before, after in events[:5]:
#             print("At sample %d stim channel went from %d to %d" % (ind, before, after))	
#    
     event_id = {'audio': 1,  'colour': 3, 'grey': 4, 'shapes': 7, 'scrambled': 8, 'button': 4096}
     color = {1: 'blue', 3: 'red', 4: 'green', 7: 'c', 8: 'black', 4096: 'yellow'}
     
     raw.save(raw_fname_new, overwrite=True)
     print events[events[:,2]==1,:].shape
     print events[events[:,2]==2,:].shape
     print events[events[:,2]==3,:].shape
     print events[events[:,2]==4,:].shape
     print events[events[:,2]==5,:].shape
     print events[events[:,2]==6,:].shape
     print events[events[:,2]==8,:].shape
     mne.viz.plot_events(events, raw.info['sfreq'], raw.first_samp, color=color, event_id=event_id)
     

#for ii, meg in enumerate(ll):
#	
#     print meg
#
#     ## Ge  t raw data
#     raw_fname1 = main_path + meg + 'block_LD1_raw.fif'
#     raw_fname = main_path + meg + 'block_LD_tsss_raw.fif'
#     raw_fname2 = main_path + meg + 'block_LD2_raw.fif'
#     raw1 = mne.io.Raw(raw_fname1, preload=True)#, preload=True
#     raw2 = mne.io.Raw(raw_fname2, preload=True)#, preload=True
#     raws=[raw1,raw2]
#     ## interpolating bad channels
#     raw1.info['bads']=bad_channels[subject_inds[ii]]
#     if raw1.info['bads']:
#         raw1.interpolate_bads(reset_bads=True)
#     print "raw1 interpolated"
#     raw2.info['bads']=bad_channels[subject_inds[ii]]
#     if raw2.info['bads']:
#         raw2.interpolate_bads(reset_bads=True)
#     print "raw2 interpolated"
#     include = []
#     exclude = []; raw1.info['bads']=[]; raw2.info['bads']=[] # bads
#     events1 = mne.find_events(raw1, stim_channel='STI101', min_duration=0.001,shortest_event=1)  
#     events2 = mne.find_events(raw2, stim_channel='STI101', min_duration=0.001,shortest_event=1) 
#     for evcnt in range(events1.shape[0]-2):
#         if events1[evcnt,2] in np.array([1,2,3,4,5]):
#             if events1[evcnt+2,2]!=16384:
#                 events1[evcnt,2]=230
#         if events1[evcnt,2] in np.array([6,7,9]):
#             if events1[evcnt+2,2]!=4096:
#                 events1[evcnt,2]=230  
#     for evcnt in range(events2.shape[0]-2):
#         if events2[evcnt,2] in np.array([1,2,3,4,5]):
#             if events2[evcnt+2,2]!=16384:
#                 events2[evcnt,2]=230
#         if events2[evcnt,2] in np.array([6,7,9]):
#             if events2[evcnt+2,2]!=4096:
#                 events2[evcnt,2]=230   
#     eventss=[events1,events2]
#     raw,events=mne.concatenate_raws(raws,events_list=eventss)
#     raw.save(raw_fname, overwrite=True)
#     mne.write_events(raw_fname[:-4]+'-eve.fif', events)
#	
#     for ind, before, after in events[:5]:
#             print("At sample %d stim channel went from %d to %d" % (ind, before, after))	
#    
#     event_id = {'visual': 1, 'hear': 2, 'hand': 3, 'neutral': 4, 'emotional': 5,'pwordc': 6, 'pworda': 7, 'pword': 9}#, 'emotional': 5,'target': 7}
#     color = {1: 'green', 2: 'yellow', 3: 'red', 4: 'c', 5: 'black', 6: 'red' ,7: 'blue', 9: 'c'}
#     mne.viz.plot_events(events, raw.info['sfreq'], raw.first_samp, color=color, event_id=event_id)

