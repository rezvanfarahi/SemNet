# -*- coding: utf-8 -*-
"""
Created on Mon May 11 12:05:11 2015

@author: rf02
"""
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
sys.path.insert(1,'/imaging/local/software/mne_python/latest_v0.9')


# for qsub
# (add ! here if needed) /imaging/local/software/anaconda/latest/x86_64/bin/python
sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
###


import numpy as np
import mne
from mne import io
from mne.minimum_norm import apply_inverse, read_inverse_operator
import matplotlib.pyplot as plt
###############################################################################
# Set parameters
data_path = '/imaging/rf02/TypLexMEG/'	# where subdirs for MEG data are
inv_path = '/imaging/rf02/TypLexMEG/'
inv_fname = 'InverseOperator_fft_1_48_clean_ica_MEG-inv.fif'
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
tmin, tmax = -0.5, 0.7
event_ids = {'cncrt_wrd': 1, 'abs_wrd': 2}
X=np.zeros((20484,35, 17,2))
for ii, meg in enumerate(ll):
    subject_no=ii#subject_inds[0]
    print subject_no
    raw_fname = data_path + meg + 'semloc_ssstf_fft_1_48_clean_ica_raw.fif' #clean_ssp_
    event_fname = event_path + meg + 'semloc_raw_ssstnew.txt'
    subjects_dir = data_path 

    tmin = -0.5
    tmax = 0.7  # Use a lower tmax to reduce multiple comparisons

	#   Setup for reading the raw data
    raw = io.Raw(raw_fname)
    events = mne.read_events(event_fname)
    stim_delay = 0.034 # delay in s
    event_id={'cncrt_wrd': 1, 'abs_wrd': 2}
    events[:,0] += np.round( raw.info['sfreq']*stim_delay )   
	###############################################################################
	# Read epochs for all channels, removing a bad one
    print "epochs"
    picks = mne.pick_types(raw.info, eeg=True, meg=True, eog=False, exclude='bads')
    if subject_no==3:
        reject = dict(eeg=150e-6, grad=200e-12, mag=4e-12) #eog=150e-6,
    else:
        reject = dict(eeg=120e-6, grad=200e-12, mag=4e-12) #eog=150e-6,
    epochs = mne.Epochs(raw, events, event_id, tmin, tmax, picks=picks, baseline=(None, 0), proj=True, reject=reject, preload=True)
    evoked1 = epochs.average()

#ed=evok.data
#em=np.abs(ed)
#emax=em.max(axis=1)
#bm=np.sqrt(np.sum(ed[:,100:500]**2,axis=1))/400
#np.mean(emax/bm)

### for epochs
    """
    epdata=epochs.get_data()
    em=np.abs(epdata)
    #emax=np.mean(em[:,:,600:900],axis=2)#
    sa=em[:,:,550:950]#.max(axis=2)
    #bm=np.sqrt(np.mean(epdata[:,:,100:500]**2,axis=2))
    bv=np.var(em[:,:,100:500],axis=2)
    edv=sa.copy()
    for ii in range(sa.shape[2]):
        edv[:,:,ii]=(sa[:,:,ii]**2)/bv
    
    ediv=np.mean(np.sqrt(edv),axis=2)
    print ediv.mean()

    #ax1 = plt.subplots(nrows=1, figsize=(6,10))
    fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(6,10))
    vmin=np.min([0, ediv.mean()-3*ediv.std()])#ediv.min()#
    vmax=np.min([ediv.max(), ediv.mean()+3*ediv.std()])#ediv.max()#
    #subplot(2,1,1)
    x1=0; x2=epdata.shape[1]; y1=0; y2=epdata.shape[0]
    im1=ax1.imshow(ediv, vmin=vmin, vmax=vmax, cmap='hot', extent=[x1,x2,y1,y2], aspect=(x2-x1)/(y2-y1))
    ax1.set_xlabel('Sensor number')
    ax1.set_ylabel('Trial number')
    
    x1=ediv.min(); x2=ediv.max(); y1=0; y2=epdata.shape[0]*epdata.shape[1]
    ax2.hist(ediv.flatten(),50)
    ax2.set_xlabel('SNR')
    ax2.set_ylabel('Sensor x Trial number')
    ax2.set_title('mean SNR:  ' + str(ediv.mean())) 
    fig.colorbar(im1, ax=ax1)
    fig.suptitle('Epochs   subject   '+meg[:10], fontsize=16)
    #plt.show(block=False)
    fig_out=data_path + meg + 'trial_SNR.jpg'
    fig.savefig(fig_out, dpi=fig.dpi)
    """
    
### for evoked
    epdata=evoked1.data
    em=np.abs(epdata)
    #emax=np.mean(em[:,:,600:900],axis=2)#
    sa=em[:,550:950]#.max(axis=2)
    #bm=np.sqrt(np.mean(epdata[:,:,100:500]**2,axis=2))
    bv=np.var(em[:,100:500],axis=1)
    edv=sa.copy()
    for ii in range(sa.shape[1]):
        edv[:,ii]=(sa[:,ii]**2)/bv
    
    ediv=np.sqrt(edv)
    print ediv.mean()
    #ax1 = plt.subplots(nrows=1, figsize=(6,10))
    fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(6,10))
    vmin=ediv.min()#np.min(0, ediv.mean()-3*ediv.std())#
    vmax=ediv.max()#ediv.mean()+3*ediv.std()#
    #subplot(2,1,1)
    x1=50; x2=50+ediv.shape[1]; y1=0; y2=ediv.shape[0]
    im1=ax1.imshow(ediv, vmin=vmin, vmax=vmax, cmap='hot', extent=[x1,x2,y1,y2], aspect=(x2-x1)/(y2-y1))
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Sensor number')
    
    x1=ediv.min(); x2=ediv.max(); y1=0; y2=ediv.shape[0]*ediv.shape[1]
    ax2.hist(ediv.flatten(),50)
    ax2.set_xlabel('SNR')
    ax2.set_ylabel('Sensor x Time Points number')
    ax2.set_title('mean SNR:  ' + str(ediv.mean())) 
    fig.colorbar(im1, ax=ax1)
    fig.suptitle('Evoked subject   '+meg[:10], fontsize=16)
    #plt.show(block=False)
    fig_out=data_path + meg + 'evoked_SNR.jpg'
    fig.savefig(fig_out, dpi=fig.dpi)
    
