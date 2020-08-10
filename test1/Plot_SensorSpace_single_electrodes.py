"""
=========================================================
Time frequency : Induced power and inter-trial phase-lock
=========================================================

This script shows how to compute induced power and inter-trial
phase-lock for a list of epochs read in a raw file given
a list of events.

"""
# Authors: Alexandre Gramfort <gramfort@nmr.mgh.harvard.edu>
#
# License: BSD (3-clause)

print __doc__

import numpy as np
import pylab as pl
import scipy.io
from mne import filter
import mne
from mne import fiff
from mne.fiff import concatenate_raws
from mne import concatenate_events, find_events
from mne.epochs import combine_event_ids
from mne.time_frequency import induced_power
from mne.datasets import sample

###############################################################################
data_path = '/imaging/nt01/MEG/Mentcalc_CueTask/Data_maxfiltered/'
data_path_out = '/imaging/nt01/MEG/Mentcalc_CueTask/Sensor_Space/Figures/'

subjects = ['meg14_0314/140724/',
            'meg14_0316/140725/',
            'meg14_0319/140725/',
            'meg14_0322/140728/',
            'meg14_0445/140731/',
            'meg14_0446/140801/',
            'meg14_0447/140804/',
            'meg14_0448/140805/',
            'meg14_0449/140805/',
            'meg14_0450/140807/',
            'meg14_0451/140808/',
            'meg14_0452/140811/',
            'meg14_0453/140811/',
            'meg14_0455/140812/',
            'meg14_0456/140814/',
            'meg14_0457/140814/',
            'meg14_0458/140815/',
            'meg14_0459/140815/',
            'meg14_0461/140818/',
            'meg14_0464/140822/',
            'meg14_0467/140826/',
            'meg14_0468/140828/',
            'meg14_0469/140901/',
            'meg14_0470/140902/']

bad_channels = (['EEG030','EEG032','EEG039','EEG043','EEG047'], #0314
                ['EEG003','EEG012','EEG043','EEG045'], #0316
                ['EEG038'], #0319
                ['EEG009'], #0322
                ['EEG002','EEG040','EEG073'], #0445
                ['EEG002','EEG008','EEG001','EEG040','EEG043'], #0446
                ['EEG004','EEG008','EEG010','EEG019','EEG020','EEG030','EEG036','EEG044','EEG058'], #0447
                ['EEG039','EEG073'], #0448
                ['EEG028','EEG045','EEG055'], #0449
                ['EEG028','EEG029','EEG039'], #0450
                ['EEG008','EEG031','EEG039'], #0451 
                ['EEG007','EEG038'], #0452
                ['EEG037','EEG069','EEG073'], #0453
                ['EEG003','EEG007','EEG015','EEG028','EEG043','EEG045','EEG046','EEG054','EEG069','EEG070'], #0455
                ['EEG004','EEG007','EEG029','EEG057'], #0456
                ['EEG006','EEG030','EEG034','EEG038','EEG044','EEG045','EEG047','EEG071'], #0457
                ['EEG002','EEG027','EEG036'], #0458
                ['EEG008','EEG028','EEG030','EEG038','EEG043','EEG045','EEG057','EEG074'], #0459
                ['EEG006','EEG028','EEG044','EEG045','EEG046','EEG070'], #0461
                ['EEG001','EEG002','EEG036','EEG055'], #0464
                ['EEG002','EEG045','EEG047','EEG073'], #0467
                ['EEG004','EEG039','EEG043','EEG050','EEG060'], #0468
                ['EEG007','EEG043'], #0469
                ['EEG044','EEG045','EEG074']) #0470

electrodes = [#'EEG010',
              'EEG011']
	      #'EEG013']
              #'EEG014',
              #'EEG016',
              #'EEG021']
              #'EEG022',
              #'EEG023',
	      #'EEG024',
	      #'EEG025',
	      #'EEG031',
	      #'EEG032',
	      #'EEG033',
	      #'EEG034',
	      #'EEG035',
	      #'EEG036',
	      #'EEG037',
	      #'EEG066',
	      #'EEG067',
	      #'EEG068']
	      #'EEG069',
	      #'EEG070']


for EEG in electrodes:
  
 power_all_cue = [];
 phase_lock_all_cue = [];
 evoked_all_cue = [];

 power_all_task = [];
 phase_lock_all_task = [];
 evoked_all_task = [];
 

 for index, meg in enumerate(subjects):

   #read concatenated events and raw file from disc 
   directory = data_path + subjects[index]
   
   raw_fname = directory + 'raw_postICA_EMEG.fif'
   #raw_fname = directory + 'session1_raw_trans1st.fif'
   
   print "reading raw file"
   raw = fiff.Raw(raw_fname, preload=True)
      
   #events_fname = directory + 'events_all-eve.fif'
   #events = mne.read_events(events_fname)
   print "Reading events from" + raw_fname
   events = mne.find_events(raw, stim_channel='STI101', min_duration=0.005) 
   #mne.write_events(directory + 'Events.txt', events)  
   
   stim_delay = 0.034 # delay in s
   events[:,0] += np.round( raw.info['sfreq']*stim_delay )

   ##################################################################################

   #extracting epochs
   tmin, tmax = -0.5, 1.0
   event_ids_cue = {'EasyMul': 101, 'ComplexMul': 102, 'EasySeq': 103, 'ComplexSeq': 104}
   event_ids_task = {'EasyMul': 131, 'ComplexMul': 132, 'EasySeq': 133, 'ComplexSeq': 134}
 
   exclude = bad_channels[index] 
   include = []
   picks = mne.pick_types(raw.info, meg=True, eeg=True, eog=True, stim=True, include=include, exclude=exclude)
   
   #print "Filtering"
   #raw.filter(l_freq=1, h_freq=40, picks=None, method='fft', l_trans_bandwidth=0.5, h_trans_bandwidth=5)
   
   print "Epoching"
   epochs_cue = mne.Epochs(raw, events, event_ids_cue, tmin, tmax, picks=picks, proj=True, baseline=(None, 0), reject=dict(grad=4000e-13, mag=4e-12, eeg=120e-6))
   epochs_task = mne.Epochs(raw, events, event_ids_task, tmin, tmax, picks=picks, baseline=(None, 0), reject=dict(grad=4000e-13, mag=4e-12, eeg=120e-6))

   #pick specific electrode
   pick_channelA = fiff.pick_channels(epochs_cue.ch_names, include=[EEG], exclude=[]);
   pick_channelB = pick_channelA + 1

    
   ##########
   #procedural tasks
   data_cue = epochs_cue.get_data()
   
   #save number of epochs
   #nr_epochs_calc = len(data_calc);
   #out_path = data_path + meg +'Num_epochs_calc.mat';
   #ep_calc = np.arange(10);
   #scipy.io.savemat(out_path, mdict={'ep_calc': nr_epochs_calc});

   
   evoked_cue = epochs_cue.average()  # compute evoked fields
   times = 1e3 * epochs_cue.times[33:467]  # change unit to ms
   evoked_cue = evoked_cue.data * 1e6  # change unit to fT / cm

   data_cue = data_cue[:, pick_channelA:pick_channelB, :]
   evoked_cue = evoked_cue[pick_channelA:pick_channelB, :]

   frequencies = np.arange(5, 35, 3)  # define frequencies of interest
   Fs = raw.info['sfreq']  # sampling in Hz
   n_cycles = frequencies / float(7)
   power_cue, phase_lock_cue = induced_power(data_cue, Fs=Fs, frequencies=frequencies,
                                 n_cycles=n_cycles, n_jobs=1, use_fft=False,
                                 zero_mean=True)

   
   #baseline corrections with ratio
   power_cue /= np.mean(power_cue[:, :, 33:167], axis=2)[:, :, None]   

   #make array with power for all subjects; average
   power_all_cue.append(power_cue[:,:,33:467])     

   #make array with phase lock for all subjects; average
   phase_lock_all_cue.append(phase_lock_cue[:,:,33:467]) 
   
   #make array with evoked data for all subjects; average
   evoked_all_cue.extend(evoked_cue[:,33:467])
 

   #########
   #retrieval tasks
   data_task= epochs_task.get_data()
   
   #save number of epochs
   #nr_epochs_retr = len(data_retr);
   #out_path = data_path + meg +'Num_epochs_retr.mat';
   #ep_retr = np.arange(10);
   #scipy.io.savemat(out_path, mdict={'ep_retr': nr_epochs_retr});

   
   evoked_task = epochs_task.average()  # compute evoked fields
   evoked_task = evoked_task.data * 1e6  # change unit to fT / cm

   data_task = data_task[:, pick_channelA:pick_channelB, :]
   evoked_task = evoked_task[pick_channelA:pick_channelB, :]

   power_task, phase_lock_task = induced_power(data_task, Fs=Fs, frequencies=frequencies,
                                 n_cycles=n_cycles, n_jobs=1, use_fft=False,
                                 zero_mean=True)

   # baseline corrections with ratio
   power_task /= np.mean(power_task[:, :, 33:167], axis=2)[:, :, None]   

   #make array with power for all subjects; average
   power_all_task.append(power_task[:,:,33:467])   

   #make array with phase lock for all subjects; average
   phase_lock_all_task.append(phase_lock_task[:,:,33:467]) 

   #make array with evoked data for all subjects; average
   evoked_all_task.extend(evoked_task[:,33:467])
   

    
   ###end loop across subjects
 
 
 #average data of all subjects
 power_average_cue = np.mean(power_all_cue, axis=0);
 phase_lock_average_cue = np.mean(phase_lock_all_cue, axis=0);
 evoked_average_cue = np.mean(evoked_all_cue, axis=0);
 mean_avg_cue = np.mean(power_average_cue[0,:,:], axis=0)
 
 power_average_task = np.mean(power_all_task, axis=0);
 phase_lock_average_task = np.mean(phase_lock_all_task, axis=0);
 evoked_average_task = np.mean(evoked_all_task, axis=0);
 mean_avg_task = np.mean(power_average_task[0,:,:], axis=0)



 ##########################################################################
 # plot averages for cue events
 import pylab as pl
 pl.clf()
 pl.subplots_adjust(0.1, 0.08, 0.96, 0.94, 0.2, 0.63)
 pl.subplot(3, 1, 1)
 pl.plot(times, mean_avg_cue)
 pl.title('Averaged Induced Power, Cue Period (%s)' % epochs_cue.ch_names[pick_channelA])
 pl.xlabel('time (ms)')
 pl.ylabel('Power')
 pl.xlim(times[0], times[-1])

 pl.subplot(3, 1, 2)
 pl.imshow(power_average_cue[0], extent=[times[0], times[-1],
                                      frequencies[0], frequencies[-1]],
          aspect='auto', origin='lower')
 pl.xlabel('Time (s)')
 pl.ylabel('Frequency (Hz)')
 pl.title('Induced power (%s)' % epochs_cue.ch_names[pick_channelA])
 pl.colorbar()

 pl.subplot(3, 1, 3)
 pl.imshow(phase_lock_average_cue[0], extent=[times[0], times[-1],
                              frequencies[0], frequencies[-1]],
          aspect='auto', origin='lower')
 pl.xlabel('Time (s)')
 pl.ylabel('Frequency (Hz)')
 pl.title('Phase-lock (%s)' % epochs_cue.ch_names[pick_channelA])
 pl.colorbar()
 #pl.show()
   
 output_figure = data_path_out + 'CuePeriod_EEGelectrode_' + epochs_cue.ch_names[pick_channelA] + '.png'
 pl.savefig(output_figure)

 
 ##############
 # plot averages for task period
 import pylab as pl
 pl.clf()
 pl.subplots_adjust(0.1, 0.08, 0.96, 0.94, 0.2, 0.63)
 pl.subplot(3, 1, 1)
 pl.plot(times, mean_avg_task)
 pl.title('Averaged Induced Power, Task Period (%s)' % epochs_task.ch_names[pick_channelA])
 pl.xlabel('time (ms)')
 pl.ylabel('Power')
 pl.xlim(times[0], times[-1])

 pl.subplot(3, 1, 2)
 pl.imshow(power_average_task[0], extent=[times[0], times[-1],
                                      frequencies[0], frequencies[-1]],
          aspect='auto', origin='lower')
 pl.xlabel('Time (s)')
 pl.ylabel('Frequency (Hz)')
 pl.title('Induced power (%s)' % epochs_task.ch_names[pick_channelA])
 pl.colorbar()

 pl.subplot(3, 1, 3)
 pl.imshow(phase_lock_average_task[0], extent=[times[0], times[-1],
                              frequencies[0], frequencies[-1]],
          aspect='auto', origin='lower')
 pl.xlabel('Time (s)')
 pl.ylabel('Frequency (Hz)')
 pl.title('Phase-lock (%s)' % epochs_task.ch_names[pick_channelA])
 pl.colorbar()
 #pl.show()
 
 output_figure = data_path_out + 'TaskPeriod_EEGelectrode_' + epochs_task.ch_names[pick_channelA] + '.png'
 pl.savefig(output_figure)
 

 ###end loop across electrode
