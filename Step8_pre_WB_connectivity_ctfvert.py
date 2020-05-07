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
import scipy
from mne.minimum_norm import  read_inverse_operator, apply_inverse_epochs, psf_ctf
from mne.connectivity import seed_target_indices, spectral_connectivity
from mne.epochs import equalize_epoch_counts
from mne import io
from scipy.stats.mstats import zscore

import os
import pickle
sys.path.insert(1,'/imaging/local/software/python_packages/scikit-learn/v0.13.1/lib.linux-x86_64-2.7/sklearn/externals')
import joblib


###############################################################################
data_path = '/imaging/rf02/Semnet/' # root directory for your MEG data
orig_path=data_path
subjects_dir = '/imaging/rf02/Semnet/'    # where your MRI subdirectories are
os.chdir('/imaging/rf02/Semnet/')
label_path='/imaging/rf02/Semnet/stc/localisers/'

# get indices for subjects to be processed from command line input
# 
print sys.argv
subject_inds = []#
for ss in sys.argv[1:]:
    subject_inds.append( int( ss ) )

subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]
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
subjects=  ['MRI_meg16_0030' ,#0
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

#labellist_path = [label_path+'col_shape_c5-lh.label',label_path+'col_shape_c5-rh.label',label_path+'auditory_c-lh.label',label_path+'auditory_c-rh.label',label_path+'hand_c-lh.label',label_path+'hand_c-rh.label']#[ 'lh.posterior_fusiform','rh.posterior_fusiform', 'lh.lateraloccipital','rh.lateraloccipital','lh.precentral','rh.precentral', 'lh.ATL_latmed','rh.ATL_latmed']#labellist = [ 'leftATL-lh','rightATL-rh']#, 'medialtempright-rh','medialtempleft-lh']'lh.postcentral','rh.postcentral',
labellist_path = [label_path+'ATL_ventromedial-lh.label',label_path+'ATL_dorsolateral-lh.label',label_path+'AG_anterior-lh.label',label_path+'AG_posterior-lh.label',label_path+'hand_postcentral-lh.label',label_path+'hand_precentral-lh.label',label_path+'col_shape_c5-lh.label']#[ 'lh.posterior_fusiform','rh.posterior_fusiform', 'lh.lateraloccipital','rh.lateraloccipital','lh.precentral','rh.precentral', 'lh.ATL_latmed','rh.ATL_latmed']#labellist = [ 'leftATL-lh','rightATL-rh']#, 'medialtempright-rh','medialtempleft-lh']'lh.postcentral','rh.postcentral',

labellist_spokes=[mne.read_label(label) for label in labellist_path]
for thisl in labellist_spokes:
    thisl.values.fill(1.0)

#fsaverage_path='/imaging/rf02/TypLexMEG'
#labelss = mne.read_labels_from_annot(subject='fsaverage', parc='rez_aparc_19oct16', subjects_dir=fsaverage_path)
##labelss.remove(labelss[-3])
##labelss.remove(labelss[-2])
##labelss.remove(labelss[-1])
#for lbl in labelss:
#    lbl.values.fill(1.0)
#label_names=[label.name for label in labelss]
#label_index=np.array([label_names.index('ATL_latmed-lh'),label_names.index('AG_hsmetafinal5-lh'),label_names.index('posterior_inferiorfrontal-lh'),label_names.index('middle_temporal-lh')])
#label_index=np.array([label_names.index('posterior_inferiorfrontal-lh'),label_names.index('middle_temporal-lh')])
#label_names.index('posterior_inferiorfrontal-lh'),label_names.index('middle_temporal-lh')
#labellist_hubs=[labelss[ii] for ii in label_index]
#labellist=labellist_hubs+labellist_spokes
labellist=labellist_spokes

print "subjects:"
print ll
n_subjects=19
stim_delay = 0.034 # delay in s ((NOTE! not in the example, taken from Olaf's script. Rezvan))
tmin, tmax = -0.3, 0.6
ctf_stc_all=range(n_subjects)
for ii, meg in enumerate(ll):
    print subject_inds[ii]
    subject_no=subject_inds[ii]
    #cov_path=data_path+meg+'noise_cov'
    #noise_cov=mne.read_cov(cov_path)
    raw_fname = data_path + meg+ 'SemDec_blocks_tsss_filt_pnt1_48_ica_raw.fif' #clean_ssp_
    event_fname = data_path + meg + 'SemDec_blocks_tsss_filt_pnt1_48_ica_raw-eve.fif'
    subjects_dir = data_path 
    
    tmin = -0.3
    tmax = 0.6  # Use a lower tmax to reduce multiple comparisons
    #   Setup for reading the raw data
    raw = io.Raw(raw_fname, preload=True)
    events = mne.read_events(event_fname)
    stim_delay = 0.034 # delay in s
    events[:,0] = events[:,0]+np.round( raw.info['sfreq']*stim_delay )
    ###############################################################################
    # Read epochs for all channels, removing a bad one
    print "epochs"
    picks = mne.pick_types(raw.info, eeg=True, meg=True, eog=False, exclude='bads')
    print "Epoching"
    if subject_inds[ii] in np.array([3,5]):#np.array([4,8,9]):
        reject = dict(eeg=150e-6, grad=200e-12, mag=4e-12)
    else:
        reject = dict(eeg=120e-6, grad=200e-12, mag=4e-12)#, eog=150e-6)
    event_ids = {'visual': 1, 'hear': 2, 'hand': 3, 'pwordc': 6}#'neutral': 4, 'emotional': 5,
    epochs = mne.Epochs(raw, events, event_ids, tmin, tmax, picks=picks, proj=True, baseline=(None, 0), reject=reject)
#    for eegcnt in range(71):
#        if eegcnt<10:
#            thiseegbad=sum(epochs.drop_log,[]).count('EEG00'+str(eegcnt))
#            if thiseegbad>=90:
#                raw.info['bads'].append(u'EEG00'+str(eegcnt))                
#        else:
#            thiseegbad=sum(epochs.drop_log,[]).count('EEG0'+str(eegcnt))
#            if thiseegbad>=90:
#                raw.info['bads'].append(u'EEG0'+str(eegcnt))
    picks = mne.pick_types(raw.info, eeg=True, meg=True, eog=False, exclude='bads')
    epochs_list=range(4)
    eid=-1
    for event_id in np.array([1,2,3,6]):
        eid = eid+1  
        epochs_list[eid]=epochs[epochs.events[:,2]==event_id]
    equalize_epoch_counts(epochs_list,method='mintime')

                   
                   
                   
	#    Equalize trial counts to eliminate bias (which would otherwise be
	#    introduced by the abs() performed below)

	###############################################################################
    print "Transform to source space"
    
    fname_inv = data_path + meg + 'InvOp_ico5oldreg_fft_pnt1_48_clean_ica_EMEG-inv.fif'
    
    method = "MNE"  # use dSPM method (could also be MNE or sLORETA)
    inverse_operator = read_inverse_operator(fname_inv)
    srcin='/imaging/rf02/TypLexMEG/fsaverage/bem/fsaverage_dist-ico-5-src.fif'
    src = mne.read_source_spaces(srcin)
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
#    snr = 3.0#np.mean(ediv)
#    print snr
    
    subject_from = subjects[subject_inds[ii]]
    subject_to = 'fsaverage'
    vertices_to = [np.arange(10242), np.arange(10242)]
#    stc_list=range(len(epochs_list))
    event_names = ['Visual', 'Hear', 'Hand', 'Pwordc']
    fname_fwd = data_path + meg + 'ico5_forward_5-3L-EMEG-fwd.fif'
    
    print "reading forward solution"
    forward_meeg = mne.read_forward_solution(fname_fwd,surf_ori=True, force_fixed=True)
    snr = 3.0#ediv.mean()
    lambda2 = 1.0 / snr ** 2
    method = 'MNE'  # can be 'MNE' or 'sLORETA'
    mode = 'svd'
    n_svd_comp = 1
    labels = [labellist[jj].morph(subject_from='fsaverage', subject_to=subjects[subject_no], subjects_dir=data_path) for jj in range(len(labellist))]
    stc_ctf_mne=psf_ctf.cross_talk_function(inverse_operator, forward_meeg, labels,method=method, lambda2=lambda2, signed=False, mode=mode, n_svd_comp=1, verbose=None)
    
    subject_from = subjects[subject_no]
    subject_to = 'fsaverage'
    vertices_to = [np.arange(10242), np.arange(10242)]
    morphmat1=mne.compute_morph_matrix(subject_from, subject_to, stc_ctf_mne.vertices, vertices_to, subjects_dir=data_path)
    ctf_stc_all[ii]=stc_ctf_mne.morph_precomputed(subject_to,vertices_to,morphmat1, subject_from)
stc_to=np.mean(ctf_stc_all)
Matx=stc_to.data[:,:-1]# Nv x Nlabels+1
Matx_zscore=np.zeros(Matx.shape)
#    for cntr in range(Matx.shape[1]-1):
#        Matx_normalised[:,cntr]=np.divide(Matx[:,cntr],Matx[:,-1])
Matx_zscore=Matx.copy()#zscore(Matx,axis=1)
assign_vertex=np.argmax(Matx_zscore, axis=1)
msort=np.sort(Matx_zscore,axis=1)
mtokeep=msort[:,-1]-msort[:,-2]
for avii in range(len(assign_vertex)):
    if Matx[avii,assign_vertex[avii]]<np.max(Matx[:,assign_vertex[avii]])/2.:
        assign_vertex[avii]=-100
        
#    assign_vertex[np.logical_or(msort[:,-1]<3, mtokeep<=1)]=-100
print "ctf finished!"

for lii,label1 in enumerate(labellist):#zip(labellist,label_index):
    print label1
    # n_cycles[np.logical_and((frequencies>10), (frequencies<=20))] = 3
    subject_from = subjects[subject_inds[ii]]
    subject_to = 'fsaverage'
    vertices_to = [np.arange(10242), np.arange(10242)]

    method = "MNE"  # use dSPM method (could also be MNE or sLORETA)

    #fname_label = label_path +  label_name + '.label'
    #label1 = mne.read_label(fname_label)
    actual_label=stc_to.in_label(label1)
    if label1.name.endswith('-lh'):
        actual_label_vert=actual_label.vertices[0]
    else:
        actual_label_vert=actual_label.vertices[1]+len(stc_to.lh_vertno)
        
    best_verts=np.intersect1d(np.where(assign_vertex==lii)[0],actual_label_vert)
    bestvert_values=Matx_zscore[best_verts,lii]
    best_101thvert=sorted(bestvert_values,reverse=True)[11]
    best_100vert=best_verts[np.where(bestvert_values>best_101thvert)]#best_verts#
    print Matx_zscore[best_100vert,lii].mean()
    label_stc_data=np.zeros((stc_to.data.shape[0],1))
    label_stc_data[best_100vert]=1.
    label_stc=mne.SourceEstimate(label_stc_data, vertices=vertices_to,tmin=1e-3 * 1., tstep=1e-3 * 1., subject='fsaverage')
    if label1.name.endswith('-lh'):
        label_new=mne.stc_to_label(label_stc, smooth=False,subjects_dir=data_path,connected=False)[0]
        label_new.name=label1.name[:-3]+'_top10verts_ctf_oldreg-lh'
        thislabelname='lh'+label1.name[:-3]
    else:
        label_new=mne.stc_to_label(label_stc, smooth=False,subjects_dir=data_path,connected=False)[1]
        label_new.name=label1.name[:-3]+'_top10verts_ctf_oldreg-rh'
        thislabelname='rh'+label1.name[:-3]
    label_new_path=label_path+str(lii)+'_'+label_new.name
    label_new.save(label_new_path)
    print label_new.vertices.shape; print "label tc extracted"
