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
subject_inds = [5]
for ss in sys.argv[1:]:
 subject_inds.append( int( ss ) )
#/imaging/rf02/Semnet/meg16_0022/160208/
#subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
print "subject_inds:"
print subject_inds
print "No rejection"

list_all = ['meg16_0022/160208/',#,
'meg16_0030/160216/',
'meg16_0032/160218/',
'meg16_0033/160218/',
'meg16_0034/160219/',
'meg16_0047/160304/'

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

     ## Ge  t raw data
     raw_fname = main_path + meg + 'block_localisers_raw_sss.fif'
     raw = mne.io.Raw(raw_fname, preload=True)#, preload=True
     ## interpolating bad channels
     if raw.info['bads']:
         raw.interpolate_bads(reset_bads=True)
     print "raw loaded"
	
     include = []
     exclude = []; raw.info['bads']=[] # bads
     picks = mne.pick_types(raw.info, meg=True, eeg=True, eog=False,
                                stim=False, include=include, exclude=exclude)
     print "picks"
     events = mne.find_events(raw, stim_channel='STI101',min_duration=0.001,shortest_event=1)    
     mne.write_events(raw_fname[:-4]+'-eve.fif', events)
	
     for ind, before, after in events[:5]:
             print("At sample %d stim channel went from %d to %d" % (ind, before, after))	
    
     event_id = {'beeps': 1, 'colour': 3, 'grey': 4, 'intact': 7, 'scrambled': 8, 'button1': 4096}
     color = {1: 'green', 3: 'yellow', 4: 'red', 7: 'c', 8: 'black', 4096: 'blue'}
     mne.viz.plot_events(events, raw.info['sfreq'], raw.first_samp, color=color, event_id=event_id)

#for ii, meg in enumerate(ll):
#	
#     print meg
#
#     ## Ge  t raw data
#     raw_fname1 = main_path + meg + 'blockld1_sss_raw.fif'
#     raw_fname = main_path + meg + 'blockld_sss_raw.fif'
#     raw_fname2 = main_path + meg + 'blockld2_sss_raw.fif'
#     raw1 = mne.io.Raw(raw_fname1, preload=True)#, preload=True
#     raw2 = mne.io.Raw(raw_fname2, preload=True)#, preload=True
#     raws=[raw1,raw2]
#     ## interpolating bad channels
#     if raw1.info['bads']:
#         raw1.interpolate_bads(reset_bads=True)
#     print "raw loaded"
#     if raw2.info['bads']:
#         raw2.interpolate_bads(reset_bads=True)
#	
#     include = []
#     exclude = []; raw1.info['bads']=[]; raw2.info['bads']=[] # bads
#     events1 = mne.find_events(raw1, stim_channel='STI101')  
#     events2 = mne.find_events(raw2, stim_channel='STI101') 
#     eventss=[events1,events2]
#          
#     raw,events=mne.concatenate_raws(raws,events_list=eventss)
#     raw.save(raw_fname, overwrite=True)
#     mne.write_events(raw_fname[:-4]+'-eve.fif', events)
#	
#     for ind, before, after in events[:5]:
#             print("At sample %d stim channel went from %d to %d" % (ind, before, after))	
#    
#     event_id = {'visual': 1, 'hear': 2, 'hand': 3, 'pword': 6}#, 'emotional': 5,'target': 7}
#     color = {1: 'green', 2: 'yellow', 3: 'red', 6: 'c'}#, 5: 'black', 7: 'blue'}
#     mne.viz.plot_events(events, raw.info['sfreq'], raw.first_samp, color=color, event_id=event_id)
