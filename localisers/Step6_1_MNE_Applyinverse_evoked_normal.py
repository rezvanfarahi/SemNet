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

print __doc__

import sys
sys.path.insert(1,'/imaging/local/software/mne_python/latest_v0.11')
sys.path.insert(1,'/imaging/rf02/scikit-learn-0.15.0')
sys.path.insert(1,'/home/rf02/.local/lib/python2.7/site-packages')
#sys.path.insert(1,'/imaging/local/software/anaconda/2.4.1/2/lib/python2.7/site-packages/scipy')
#import scipy
#print scipy.__version__
# for qsub
# (add ! here if needed) /imaging/local/software/anaconda/latest/x86_64/bin/python
#sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
###

import numpy as np
import mne

from mne.minimum_norm import apply_inverse, read_inverse_operator
from mne.epochs import equalize_epoch_counts
from mne import io
import os


###############################################################################
data_path = '/imaging/rf02/Semnet/' # root directory for your MEG data
orig_path=data_path
subjects_dir = '/imaging/rf02/Semnet/'    # where your MRI subdirectories are
os.chdir('/imaging/rf02/Semnet/')

# get indices for subjects to be processed from command line input
# 
print sys.argv
subject_inds = []#
for ss in sys.argv[1:]:
    subject_inds.append( int( ss ) )

subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13]
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


# subjects names used for MRI data
subjects=  [
            'MRI_meg16_0045' ,#8
            'MRI_meg16_0052' ,#10
            'MRI_meg16_0056' ,#11
    	       'MRI_meg16_0069' ,#12
 	       'MRI_meg16_0070' ,#13
            'MRI_meg16_0072' ,#15
            'MRI_meg16_0073' ,#16
            'MRI_meg16_0075' ,#17
            'MRI_meg16_0078' ,#18
            'MRI_meg16_0082' ,#19
            'MRI_meg16_0086' ,#20
            'MRI_meg16_0097' ,#21
            'MRI_meg16_0122' ,#22
            'MRI_meg16_0125' ,#24
            ]

ll = []
for ss in subject_inds:
 ll.append(list_all[ss])
event_names = ['audio',  'colour', 'grey', 'shapes', 'scrambled']
stc_alllists=[list() for lii in range(len(event_names))]
print "subjects:"
print ll
n_subjects=14
stim_delay = 0.034 # delay in s ((NOTE! not in the example, taken from Olaf's script. Rezvan))
tmin, tmax = -0.5, 0.3
for ii, meg in enumerate(ll):
    print subject_inds[ii]
    #cov_path=data_path+meg+'noise_cov'
    #noise_cov=mne.read_cov(cov_path)
    raw_fname = data_path + meg+ 'block_localisers_tsss_filt_pnt1_48_ica_raw.fif' #clean_ssp_
    event_fname = orig_path + meg + 'block_localisers_tsss_ds_raw-eve.fif'
    subjects_dir = data_path 
    
    tmin = -0.1
    tmax = 0.5  # Use a lower tmax to reduce multiple comparisons
    #   Setup for reading the raw data
    raw = io.Raw(raw_fname)
    events = mne.read_events(event_fname)
    stim_delay = 0.034 # delay in s
    events[:,0] = events[:,0]+np.round(raw.info['sfreq']*stim_delay )#raw.info['sfreq'] instead of 1000?
    ###############################################################################
    # Read epochs for all channels, removing a bad one
    print "epochs"
    picks = mne.pick_types(raw.info, eeg=True, meg=True, eog=False, exclude='bads')
    print "Epoching"
    if subject_inds[ii] in np.array([0]):#np.array([4,8,9]):
        reject = dict(eeg=150e-6, grad=200e-12, mag=4e-12)
    else:
        reject = dict(eeg=120e-6, grad=200e-12, mag=4e-12)#, eog=150e-6)
    event_ids = {'audio': 1,  'colour': 3, 'grey': 4, 'shapes': 7, 'scrambled': 8}#4096# for button
    epochs = mne.Epochs(raw, events, event_ids, tmin, tmax, picks=picks, proj=True, baseline=(None,0), reject=reject) #-0.43,-0.33 for button
    for eegcnt in range(71):
        if eegcnt<10:
            thiseegbad=sum(epochs.drop_log,[]).count('EEG00'+str(eegcnt))
            if thiseegbad>=90:
                raw.info['bads'].append(u'EEG00'+str(eegcnt))                
        else:
            thiseegbad=sum(epochs.drop_log,[]).count('EEG0'+str(eegcnt))
            if thiseegbad>=90:
                raw.info['bads'].append(u'EEG0'+str(eegcnt))
    picks = mne.pick_types(raw.info, eeg=True, meg=True, eog=False, exclude='bads')
    epochs_list=range(5)
    eid=-1
    for event_id in np.array([1,3,4,7,8]):#np.array([4096]):#
        eid=eid+1
        epochs_list[eid]=epochs[epochs.events[:,2]==event_id]
#    equalize_epoch_counts(epochs_list,method='mintime')

                   
                   
                   
	#    Equalize trial counts to eliminate bias (which would otherwise be
	#    introduced by the abs() performed below)

	###############################################################################
    print "Transform to source space"
    
    fname_inv = data_path + meg + 'InvOp_ico5newreg_fft_pnt1_48_localisers_ica_EMEG-inv.fif'
    
    method = "MNE"  # use dSPM method (could also be MNE or sLORETA)
    inverse_operator = read_inverse_operator(fname_inv)
#    sample_vertices = [s['vertno'] for s in inverse_operator['src']]
    #    Let's average and compute inverse, resampling to speed things up
#    evoked=epochs.average()    
#    epdata=evoked.data 
#    em=np.abs(epdata)
#    #emax=np.mean(em[:,:,600:900],axis=2)#
#    sa=em[:,550:950]#.max(axis=2)
#    #bm=np.sqrt(np.mean(epdata[:,:,100:500]**2,axis=2))
#    bv=np.var(em[:,100:500],axis=1)
#    edv=sa.copy()
#    for ii in range(sa.shape[1]):
#        edv[:,ii]=(sa[:,ii]**2)/bv
#
#    ediv=np.sqrt(edv)
    snr = 3.0#np.mean(ediv)
    print snr
    lambda2 = 1.0 / snr ** 2
    subject_from = subjects[subject_inds[ii]]
    subject_to = 'fsaverage'
    vertices_to = [np.arange(10242), np.arange(10242)]
    vertices_avg = [np.arange(10242), np.arange(10242)]      

#    stc_list=range(len(epochs_list))
    event_names = ['audio',  'colour', 'grey', 'shapes', 'scrambled']#['button']#
    srcin='/imaging/rf02/TypLexMEG/fsaverage/bem/fsaverage_dist-ico-5-src.fif'
    src_avg = mne.read_source_spaces(srcin)
    tminl=0
    tstepl=1
    thislbl_data=np.ones((len(vertices_avg[0])+len(vertices_avg[1]),1))
    thisstc_tolabel = mne.SourceEstimate(thislbl_data, vertices=vertices_avg,tmin=1e-3 * tminl, tstep=1e-3 * tstepl, subject='fsaverage')
    thisstc_label=mne.stc_to_label(thisstc_tolabel, src_avg)
    thisscale=-np.hstack((mne.label.label_sign_flip(thisstc_label[0],src_avg),mne.label.label_sign_flip(thisstc_label[1],src_avg)))
    
#    evoked_list=range(len(epochs_list))
#    condition_list=range(len(epochs_list))
    for evcnt in range(len(epochs_list)):
        this_evoked = epochs_list[evcnt].average()
        this_condition = apply_inverse(this_evoked, inverse_operator, lambda2, method, pick_ori="normal")
        this_morphmat=mne.compute_morph_matrix(subject_from, subject_to, this_condition.vertices, vertices_to, subjects_dir=data_path)
        this_stc=this_condition.morph_precomputed(subject_to,vertices_to,this_morphmat, subject_from)
        stc_alllists[evcnt].append(this_stc)
stc_allmean=[np.mean(tstc) for tstc in stc_alllists]
tstepds=1/raw.info['sfreq']
out_path='/imaging/rf02/Semnet/stc/GrandAverage/evoked/localisers/signed/'
for evcnt, tstc in enumerate(stc_allmean):
    predata1=tstc.data.copy()
    ntimes=len(tstc.times)
    #    label_data1=np.tile(thisscale,[ntimes,1]).transpose()*label_predata
    signeddata=np.tile(thisscale[:,np.newaxis],(1,predata1.shape[1]))*predata1
    stc_signeddata = mne.SourceEstimate(signeddata, vertices=vertices_avg,tmin=1e-3 * tmin, tstep=tstepds, subject='fsaverage')
    data_out = out_path +'Grandavg_firstMorphed_ico_newreg_localisers_ica_'+event_names[evcnt]+'_Source_Evoked_m100_500' #not really -100_500 for the button
    stc_signeddata.save(data_out)
    # Now we can visualize the coherence using the plot method.
#    brain = coh_stc.plot('sample', 'inflated', 'both',
#                     time_label='Coherence %0.1f Hz',
#                     subjects_dir=subjects_dir,
#                     clim=dict(kind='value', lims=(0.25, 0.4, 0.65)))
#    brain.show_view('lateral')
#
#thisscale=mne.label.label_sign_flip(tlab,src_sub)#np.hstack((mne.label.label_sign_flip(tlab,src_sub),mne.label.label_sign_flip(thisstc_label[1],src_sub)))[winner_verts]#mne.label.label_sign_flip(this_label,src_avg)[winner_verts]#1#Matx[winner_verts,this_labelc]/np.max(Matx[winner_verts,this_labelc])
#        label_predata1=condition1.in_label(tlab).data
#        #            label_data1=np.tile(thisscale,[ntimes,1]).transpose()*label_predata
#        label_data1=np.tile(thisscale[:,np.newaxis],(1,label_predata1.shape[1]))*label_predata1
#        ini_tc1=np.mean(label_data1, axis=0)
#        label_predata2=condition2.in_label(tlab).data
#        #            label_data1=np.tile(thisscale,[ntimes,1]).transpose()*label_predata
#        label_data2=np.tile(thisscale[:,np.newaxis],(1,label_predata2.shape[1]))*label_predata2
#        ini_tc2=np.mean(label_data2, axis=0)
#        conabs_mat[tlc,:,0]=mne.extract_label_time_course(condition1,tlab,src_sub,mode='mean_flip')
#        conabs_mat[tlc,:,1]=mne.extract_label_time_course(condition2,tlab,src_sub,mode='mean_flip')
##    X[ii,:,:,:]=conabs_mat
#    out_file=out_path + meg[:10] + '_SemLoc_Evoked_5ROIs_ctfdip_exttc_'+str(np.abs(roivert))+'verts_avg.mat'
#    scio.savemat(out_file,{'conabs_mat':conabs_mat})
#
#        
#vertices_avg = [np.arange(10242), np.arange(10242)]
#    vertices_sub=[np.arange(forward['src'][0]['nuse']),np.arange(forward['src'][1]['nuse'])]
##    conind=0
##    for event_no,event_name in enumerate(['Concrete','Abstract']):
###        conind=conind+1
###        print conind
##                # Load data
##        
##        fname = data_path + meg + 'firstMorphed_ico_signed_SemLoc_ica_'+event_name+'_Source_Evoked_m500_700' 
##        this_stc = mne.read_source_estimate(fname)  
#        tmin1=-500
#        tstep1=1
#        thisstc_data=np.ones((this_stc.data.shape[0],1))
#        thisstc_tolabel = mne.SourceEstimate(thisstc_data, vertices=vertices_avg,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
#        thisstc_label=mne.stc_to_label(thisstc_tolabel, src_avg)

    #evoked1 = whiten_evoked(evoked11, noise_cov)
    #evoked1.resample(50)