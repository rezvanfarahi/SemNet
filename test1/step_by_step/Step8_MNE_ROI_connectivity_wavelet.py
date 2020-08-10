"""
=========================================================================
Compute source space connectivity and visualize it using a circular graph
=========================================================================

This example computes the all-to-all connectivity between 68 regions in
source space based on dSPM inverse solutions and a FreeSurfer cortical
parcellation. The connectivity is visualized using a circular graph which
is ordered based on the locations of the regions.
"""

# Authors: Martin Luessi <mluessi@nmr.mgh.harvard.edu>
#          Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
#          Nicolas P. Rougier (graph code borrowed from his matplotlib gallery)
#
# License: BSD (3-clause)

print(__doc__)
import sys

#sys.path.append('/imaging/local/software/python_packages/scikit-learn/v0.14.1')
sys.path.append('/imaging/local/software/mne_python/latest_v0.9')

# for qsub
# (add ! here if needed) 
sys.path.append('/imaging/local/software/EPD/latest/x86_64/bin/python')
#sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
###

import numpy as np
import mne
from mne.io import Raw
from mne.minimum_norm import (read_inverse_operator, apply_inverse, apply_inverse_epochs)
import scipy.io
from mne import find_events
from mne.connectivity import (spectral_connectivity, phase_slope_index)
import numpy as np
from mne import read_labels_from_annot


data_path = '/imaging/rf02/TypLexMEG/' # root directory for your MEG data
subjects_dir = '/imaging/rf02/TypLexMEG/'    # where your MRI subdirectories are

event_path = event_path = '/group/erp/data/olaf.hauk/MEG/TypLexM/data/MF22_new/'   # where event files are

label_path = '/imaging/rf02/TypLexMEG/fsaverage/label/'

inv_path = '/imaging/rf02/TypLexMEG/'
inv_fname = 'InverseOperator_fft_1_48_clean_ica_EMEG-inv.fif'

# get indices for subjects to be processed from command line input
# 
print sys.argv
subject_inds = []
for ss in sys.argv[1:]:
   subject_inds.append( int( ss ) )



#subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]#0,1,2,3,4,5,6,7,8]#,
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
'meg11_0147/110603/', 


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
Dmat_cncrt=np.zeros((72,72,4,17))
Dmat_abs=np.zeros((72,72,4,17))
cnt1=subjects[subject_inds[0]]
for meg in ll:

#    print meg
#    #con_cncrt=np.divide(con_cncrt,con_blcn)
#    path_c=data_path+ meg+'plv_concrete_cwtmorlet.npy'
#    plv_concrete_cwtmorlet=np.load(path_c)
#    scipy.io.savemat(path_c[:-4],{'plv_concrete_cwtmorlet':plv_concrete_cwtmorlet})
#    print "con finished"
#    path_a=data_path+ meg+'plv_abstract_cwtmorlet.npy'
#    plv_abstract_cwtmorlet=np.load(path_a)
#    scipy.io.savemat(path_a[:-4],{'plv_abstract_cwtmorlet':plv_abstract_cwtmorlet})
#    print "abs finished"
#    path_s=data_path+ meg+'plv_subtract_cwtmorlet.npy'
#    plv_subtract_cwtmorlet=np.load(path_s)
#    scipy.io.savemat(path_s[:-4],{'plv_subtract_cwtmorlet':plv_subtract_cwtmorlet})
#    print "sub finished"
#    #con_abs=np.divide(con_abs,con_blab)

    #cnt1=cnt1+1
    
    print cnt1
    print meg
    fname_inv = inv_path + meg + 'InverseOperator_fft_1_48_clean_ica_EMEG-inv.fif'
    fname_raw = data_path + meg + '/semloc_ssstf_fft_1_48_clean_ica_raw.fif'
    fname_event = event_path + meg + 'semloc_raw_ssstnew.txt'

    # Load data
    inverse_operator = read_inverse_operator(fname_inv)
    srcin='/imaging/rf02/TypLexMEG/fsaverage/bem/fsaverage-ico-5-src.fif'
    src = mne.read_source_spaces(srcin)#inverse_operator['src']
#    raw = Raw(fname_raw, preload=True)
#    
#    
#    events = mne.read_events(fname_event)
#    stim_delay = 0.034
#    events[:,0] += np.round( raw.info['sfreq']*stim_delay )
#    # Add a bad channel
#    #raw.info['bads'] += ['MEG 2443']
#
#    # Pick MEG channels
#    picks = mne.pick_types(raw.info, meg=True, eeg=True, stim=False, eog=True, exclude='bads')
#
#    # Define epochs for left-auditory condition
#    tmin, tmax = -0.5, 0.7
#    if subject_inds[0]==3:
#        reject = dict(eeg=150e-6, grad=200e-12, mag=4e-12) #eog=150e-6,
#    else:
#        reject = dict(eeg=120e-6, grad=200e-12, mag=4e-12) #eog=150e-6, 
#    event_ids = {'cncrt_wrd': 1, 'abs_wrd': 2}
#    method = "MNE"  # use dSPM method (could also be MNE or sLORETA)
#
#    print "Read epochs"
#    epochs = mne.Epochs(raw, events, event_ids, tmin, tmax, picks=picks, baseline=(None, 0), proj=True, reject=reject)
#    epochs_concrete = epochs['cncrt_wrd']
#    epochs_abstract = epochs['abs_wrd']
#    #epochs=epochs_abstract
#    events_no_concrete=epochs_concrete.copy().get_data().shape[0]
#    events_no_abstract=epochs_abstract.copy().get_data().shape[0]
#    stcs_concrete=range(events_no_concrete)
##mne.epochs.equalize_epoch_counts([epochs_concrete,epochs_abstract])
#    for ccc in range(len(stcs_concrete)):
#        print ccc		
#        data_in = data_path + meg + 'Morphed_SemLoc_epochs' + str(ccc) + '_icaclean_Concrete_Source_Evoked_m500_700'
#        stcbefore=mne.read_source_estimate(data_in)
#        stcbefore.resample(200)
#        stcs_concrete[ccc]=stcbefore
#        
#    stcs_abstract=range(events_no_abstract)
##mne.epochs.equalize_epoch_counts([epochs_concrete,epochs_abstract])
#    for ccc in range(len(stcs_abstract)):
#        print ccc		
#        data_in = data_path + meg + 'Morphed_SemLoc_epochs' + str(ccc) + '_icaclean_Abstract_Source_Evoked_m500_700'
#        stcbefore=mne.read_source_estimate(data_in)
#        stcbefore.resample(200)
#        stcs_abstract[ccc]=stcbefore
#        
#    raw.resample(200)
#	# Compute inverse solution and for each epoch. By using "return_generator=True"
#	# stcs will be a generator object instead of a list.
##	snr = 1.0  # use lower SNR for single epochs
##	lambda2 = 1.0 / snr ** 2
#    method = "MNE"  # use dSPM method (could also be MNE or sLORETA)
#    # Get labels for FreeSurfer 'aparc' cortical parcellation with 34 labels/hemi
#    labels = mne.read_labels_from_annot(subject='fsaverage', parc='rezvan_parc_mne_ctf_mn', subjects_dir=subjects_dir)
#    label_colors = [label.color for label in labels]
# 
#    ##	epochs = mne.Epochs(raw, events, event_ids, tmin, tmax, picks=picks,
#    #	#	            baseline=(None, 0), reject=reject)
#    #	#epochs_concrete = epochs['cncrt_wrd']
#    #	#epochs_abstract = epochs['abs_wrd']
#    ##mne.epochs.equalize_epoch_counts([epochs_concrete,epochs_abstract])
#    ##	stcs_concrete = apply_inverse_epochs(epochs_concrete, inverse_operator, lambda2, method, pick_ori="normal", return_generator=True)
#    #	#stcs_abstract = apply_inverse_epochs(epochs_abstract, inverse_operator, lambda2, method, pick_ori="normal", return_generator=True)
#    #	#	
#    #	# Average the source estimates within each label using sign-flips to reduce
#    #	# signal cancellations, also here we return a generator
#    
 
    #label_mne_ctf_mn_ts_concrete_resampled = mne.extract_label_time_course(stcs_concrete, labels[:-2], src, mode='pca_flip',return_generator=False)
    #label_mne_ctf_mn_ts_abstract_resampled = mne.extract_label_time_course(stcs_abstract, labels[:-2], src, mode='pca_flip',return_generator=False)
    path_c=data_path+ meg+'label_mne_ctf_mn_ts_concrete_resampled.mat'
    ltc=scipy.io.loadmat(path_c,mat_dtype=True); label_mne_ctf_mn_ts_concrete_resampled=ltc['label_mne_ctf_mn_ts_concrete_resampled']
    #scipy.io.savemat(path_c,{'label_mne_ctf_mn_ts_concrete_resampled':label_mne_ctf_mn_ts_concrete_resampled})
    path_a=data_path+ meg+'label_mne_ctf_mn_ts_abstract_resampled.mat'
    #scipy.io.savemat(path_a,{'label_mne_ctf_mn_ts_abstract_resampled':label_mne_ctf_mn_ts_abstract_resampled})
    lta=scipy.io.loadmat(path_a,mat_dtype=True); label_mne_ctf_mn_ts_abstract_resampled=lta['label_mne_ctf_mn_ts_abstract_resampled']
    print "lbl ts"

#	
  
    		
    fmin = (4., 8., 13., 31.)
    fmax = (7., 12., 30, 45.)
    tmin2 = -0.2
    tmax2 = 0.55
    tminb2 = -0.4
    tmaxb2 = 0.
    sfreq = 200.0#raw.info['sfreq']  # the sampling frequency
    con_methods = ['coh']
    #cwt_n_cycles=np.array([2,2,2,2,3,3,3,3,3,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7])
    cwt_frequencies = np.arange(4, 46, 1)#np.array([6,10,22,38])#
    cwt_n_cycles = cwt_frequencies / float(2)#np.array([2,3,5,7])#
    #cwt_n_cycles[cwt_frequencies<=20] = 2
    #cwt_n_cycles[cwt_frequencies<=20] += np.arange(0,1,0.059)#0.117)#0.059)
    #con_blcn, freqs_blcn, times_blcn, n_epochs_blcn, n_tapers_blcn = spectral_connectivity(label_ts_concrete, method=con_methods, mode='multitaper', sfreq=sfreq, fmin=fmin, fmax=fmax, faverage=True, tmin=tminb2, tmax=tmaxb2, mt_adaptive=True, n_jobs=2)
    
    #con_blab, freqs_bl, times_bl, n_epochs_bl, n_tapers_bl = spectral_connectivity(label_ts_abstract, method=con_methods, mode='multitaper', sfreq=sfreq, fmin=fmin, fmax=fmax, faverage=True, tmin=tminb2, tmax=tmaxb2, mt_adaptive=True, n_jobs=2)
    
    con_cncrt, freqs_cncrt, times_cncrt, nnn, ppp= spectral_connectivity(label_mne_ctf_mn_ts_concrete_resampled, method=con_methods, mode='cwt_morlet', cwt_n_cycles=cwt_n_cycles, cwt_frequencies=cwt_frequencies, sfreq=sfreq, fmin=fmin, fmax=fmax, faverage=True)#, n_jobs=4)
    ##con_cncrt, freqs_cncrt, times_cncrt, nnn, ppp= phase_slope_index(label_mne_ctf_mn_ts_concrete_resampled,  mode='cwt_morlet', cwt_n_cycles=cwt_n_cycles, cwt_frequencies=cwt_frequencies, sfreq=sfreq, fmin=fmin, fmax=fmax)#, faverage=True)#, n_jobs=4)

    #con_cncrt=np.divide(con_cncrt,con_blcn)
    path_c=data_path+ meg+'pca_coh_concrete_cwtmorlet_resmpled'
    np.save(path_c,con_cncrt)
    scipy.io.savemat(path_c,{'con_cncrt':con_cncrt})
    print "con finished"
    con_abs, freqs_abs, times_abs, nnn, ppp = spectral_connectivity(label_mne_ctf_mn_ts_abstract_resampled,method=con_methods, mode='cwt_morlet', cwt_n_cycles=cwt_n_cycles, cwt_frequencies=cwt_frequencies, sfreq=sfreq, fmin=fmin, fmax=fmax, faverage=True)#,  n_jobs=4)
    ##con_abs, freqs_abs, times_abs, nnn, ppp = phase_slope_index(label_mne_ctf_mn_ts_abstract, mode='cwt_morlet', cwt_n_cycles=cwt_n_cycles, cwt_frequencies=cwt_frequencies, sfreq=sfreq, fmin=fmin, fmax=fmax)#, faverage=True)#,  n_jobs=4)

    path_a=data_path+ meg+'pca_coh_abstract_cwtmorlet_resmpled'
    np.save(path_a,con_abs)
    scipy.io.savemat(path_a,{'con_abs':con_abs})
    print "abs finished"
    con_sub=con_cncrt-con_abs
    path_s=data_path+ meg+'pca_coh_subtract_cwtmorlet_resmpled'
    np.save(path_s,con_sub)
    scipy.io.savemat(path_s,{'con_sub':con_sub})
    print "sub finished"
    #con_abs=np.divide(con_abs,con_blab)
	
#"""
#	# con is a 3D array, get the connectivity for the first (and only) freq. band
#	# for each method
#	con_res = con_cncrt
#	#for method, c in zip(con_methods, con):
#	 #   con_res[method] = c[:, :, 0]
#
#	# Now, we visualize the connectivity using a circular graph layout
#
#	# First, we reorder the labels based on their location in the left hemi
#	label_names = [label.name for label in labels]
#
#	lh_labels = [name for name in label_names if name.endswith('lh')]
#
#	# Get the y-location of the label
#	label_ypos = list()
#	for name in lh_labels:
#	    idx = label_names.index(name)
#	    ypos = np.mean(labels[idx].pos[:, 1])
#	    label_ypos.append(ypos)
#
#	# Reorder the labels based on their location
#	lh_labels = [label for (ypos, label) in sorted(zip(label_ypos, lh_labels))]
#
#	# For the right hemi
#	rh_labels = [label[:-2] + 'rh' for label in lh_labels]
#
#	# Save the plot order and create a circular layout
#	node_order = list()
#	node_order.extend(lh_labels[::-1])  # reverse the order
#	node_order.extend(rh_labels)
#
#	node_angles = circular_layout(label_names, node_order, start_pos=90,
#		                      group_boundaries=[0, len(label_names) / 2])
#	Dmat_cncrt[:,:,:,cnt1]=con_cncrt
#	Dmat_abs[:,:,:,cnt1]=con_abs
#
#Dmatc=Dmat_cncrt.mean(axis=3)
#Dmata=Dmat_abs.mean(axis=3)
#
#	# Plot the graph using node colors from the FreeSurfer parcellation. We only
#	# show the 300 strongest connections.
##plot_connectivity_circle(con_res[:,:,0], label_names, n_lines=300,node_angles=node_angles, node_colors=label_colors,title='All-to-All Connectivity Concrete Word' 'Condition (coh)')
#
#import matplotlib.pyplot as plt
#	#plt.savefig('circle.png', facecolor='black')
#
#	# Plot connectivity for both methods in the same plot
#fig = plt.figure(num=None, figsize=(12, 12), facecolor='black')
##no_names = [''] * len(label_names)
#for ii in range(1):#Dmatc.shape[2]):
#	plot_connectivity_circle(Dmatc[:,:,ii], label_names, n_lines=300,node_angles=node_angles, node_colors=label_colors, title='coh', padding=0, fontsize_colorbar=4, fig=fig, subplot=(2, 2, ii + 1))
#
#fig = plt.figure(num=None, figsize=(12, 12), facecolor='black')
##no_names = [''] * len(label_names)
#for ii in range(1):#Dmata.shape[2]):
#	plot_connectivity_circle(Dmata[:,:,ii], label_names, n_lines=300,node_angles=node_angles, node_colors=label_colors, title='coh', padding=0, fontsize_colorbar=4, fig=fig, subplot=(2, 2, ii + 1))
#
#Dmats=Dmatc-Dmata;
#fig = plt.figure(num=None, figsize=(12, 12), facecolor='black')
##no_names = [''] * len(label_names)
#for ii in range(1):#Dmats.shape[2]):
#	plot_connectivity_circle(Dmats[:,:,ii], label_names, n_lines=30,node_angles=node_angles, node_colors=label_colors, title='coh', padding=0, fontsize_colorbar=8, fig=fig, subplot=(1, 1, ii + 1),colormap='jet')
#
#plt.show()
#"""
