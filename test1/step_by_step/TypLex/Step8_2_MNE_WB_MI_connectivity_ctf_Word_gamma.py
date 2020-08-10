# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 14:55:11 2015

@author: rf02
"""

# Author: Martin Luessi <mluessi@nmr.mgh.harvard.edu>
#
# License: BSD (3-clause)

print(__doc__)
import sys

#sys.path.append('/imaging/local/software/python_packages/nibabel/1.3.0')
#sys.path.append('/imaging/local/software/python_packages/pysurfer/v0.4')
#sys.path.append('/imaging/local/software/python_packages/scikit-learn/v0.14.1')
#sys.path.append('/imaging/local/software/mne_python/latest')

# for qsub
# (add ! here if needed)
sys.path.append('/imaging/local/software/EPD/latest/x86_64/bin/python')
sys.path.insert(1,'/imaging/rf02/TypLexMEG/mne_python0.9')
#sys.path.insert(1,'/imaging/local/software/mne_python/latest_v0.9full')
#/imaging/local/software/mne_python/latest
#sys.path.append('/imaging/local/software/python_packages/scikit-learn/v0.16.1')
sys.path.insert(1,'/imaging/local/software/python_packages/scikit-learn/v0.13.1/lib.linux-x86_64-2.7/sklearn/externals')
import joblib
###


import os
import numpy as np
import mne
import scipy
from scipy.stats.mstats import zscore
#import sklearn
from mne.io import Raw
from mne.minimum_norm import psf_ctf_28Oct15
from mne.minimum_norm import (read_inverse_operator,  apply_inverse, apply_inverse_epochs)
from mne.connectivity import seed_target_indices, spectral_connectivity
from numpy import *
import numpy
import time

def entropy(counts):
    '''Compute entropy.'''
    ps = counts/float(sum(counts))  # coerce to float and normalize
    ps = ps[nonzero(ps)]            # toss out zeros
    H = -sum(ps * numpy.log2(ps))   # compute entropy
    
    return H

def mi(x, y, bins):
    '''Compute mutual information'''
    counts_xy = histogram2d(x, y, bins=bins)[0]
    counts_x  = histogram(  x,    bins=bins)[0]
    counts_y  = histogram(  y,    bins=bins)[0]
    
    H_xy = entropy(counts_xy)
    H_x  = entropy(counts_x)
    H_y  = entropy(counts_y)
    
    return H_x + H_y - H_xy
data_path = '/imaging/rf02/TypLexMEG/' # root directory for your MEG data
os.chdir(data_path)
event_path = '/group/erp/data/olaf.hauk/MEG/TypLexM/data/MF22_new/'   # where event files are
label_path = '/imaging/rf02/TypLexMEG/fsaverage/label/'
inv_path = '/imaging/rf02/TypLexMEG/'
subjects_dir = data_path
print sys.argv
subject_inds = []
for ss in sys.argv[1:]:
   subject_inds.append( int( ss ) )

#subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
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




#label_name_lh = ['atlleft-lh', 'atlright-rh']


# Labels/ROIs'lh.lateralorbitofrontal','rh.lateralorbitofrontal']#
#labellist = ['lh.ATL_latmed','lh.supramarginal','lh.posterior_inferiorfrontal','lh.middle_temporal']#[ 'lh.posterior_fusiform','rh.posterior_fusiform', 'lh.lateraloccipital','rh.lateraloccipital','lh.precentral','rh.precentral', 'lh.ATL_latmed','rh.ATL_latmed']#labellist = [ 'leftATL-lh','rightATL-rh']#, 'medialtempright-rh','medialtempleft-lh']'lh.postcentral','rh.postcentral',
labelss = mne.read_labels_from_annot(subject='fsaverage', parc='rez_aparc_19oct16', subjects_dir=data_path)
#labelss.remove(labelss[-3])
#labelss.remove(labelss[-2])
for lbl in labelss:
    lbl.values.fill(1.0)
label_names=[label.name for label in labelss]
label_index=np.array([label_names.index('ATL_latmed-lh'),label_names.index('ATL_latmed-rh'),label_names.index('AG_hsmetafinal5-lh'),label_names.index('posterior_inferiorfrontal-lh'),label_names.index('middle_temporal-lh')])

labellist=[labelss[ii] for ii in label_index]#[labelss[ii] for ii in label_index]
print labellist
#mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/lh.supramarginal.label', supramarginal)
#
#posterior_inferiorfrontal=labelss[label_name.index('posterior_inferiorfrontal-lh')]
#mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/lh.posterior_inferiorfrontal.label', posterior_inferiorfrontal)
#
#middle_temporal=labelss[label_name.index('middle_temporal-lh')]
#mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/lh.middle_temporal.label', middle_temporal)

ii=-1
for meg in ll:
    print subject_inds
    subject_no=subject_inds[0]
    tmin, tmax = -0.5, 0.7
    fname_fwd = data_path + meg + 'ico5_forward_5-3L-EMEG-fwd.fif'
    forward = mne.read_forward_solution(fname_fwd,surf_ori=True, force_fixed=True)
    if subject_inds[0] in np.array([3,4,13]):
        reject = dict(eeg=150e-6, grad=200e-12, mag=4e-12) #eog=150e-6,
    else:
        reject = dict(eeg=120e-6, grad=200e-12, mag=4e-12) #eog=150e-6,
# frequency bands with variable number of cycles for wavelets
    stim_delay = 0.034 # delay in s
#  n_cycles[np.logical_and((frequencies>10), (frequencies<=20))] = 3

    method = "MNE"  # use dSPM method (could also be MNE or sLORETA)
    ii=ii+1
    fname_inv = inv_path + meg + 'typlex_InvOp_ico5diagreg_fft_1_48_clean_ica_EMEG-inv.fif'
    fname_raw = data_path + meg + 'typlex_pw_ssstf_fft_ica_gamma_raw.fif'
    fname_event = event_path + meg + 'typlex_pw_raw_ssstnew.txt'
# Load data
    inverse_operator = read_inverse_operator(fname_inv)
    srcin='/imaging/rf02/TypLexMEG/fsaverage/bem/fsaverage-ico-5-src.fif'
    src = mne.read_source_spaces(srcin)#inverse_operator['src']  # the source space used
    raw = Raw(fname_raw)
    events = mne.read_events(fname_event)
    stim_delay = 0.034 # delay in s
    events[:,0] += np.round( raw.info['sfreq']*stim_delay ) 
    
    ### select 100 most sensitive labels
    snr = 3.0#ediv.mean()
    lambda2 = 1.0 / snr ** 2
    method = 'MNE'  # can be 'MNE' or 'sLORETA'
    mode = 'svd'
    n_svd_comp = 1
    labels = [labelss[jj].morph(subject_from='fsaverage', subject_to=subjects[subject_no], subjects_dir=data_path) for jj in range(len(labelss))]
    stc_ctf_mne=psf_ctf_28Oct15.cross_talk_function(inverse_operator, forward, labels,method=method, lambda2=lambda2, signed=False, mode=mode, n_svd_comp=1, verbose=None)
    
    subject_from = subjects[subject_no]
    subject_to = 'fsaverage'
    vertices_to = [np.arange(10242), np.arange(10242)]
    morphmat1=mne.compute_morph_matrix(subject_from, subject_to, stc_ctf_mne.vertices, vertices_to, subjects_dir=data_path)
    stc_to=stc_ctf_mne.morph_precomputed(subject_to,vertices_to,morphmat1, subject_from)
    Matx=stc_to.data# Nv x Nlabels+1
    Matx_zscore=np.zeros(Matx[:,:72].shape)
#    for cntr in range(Matx.shape[1]-1):
#        Matx_normalised[:,cntr]=np.divide(Matx[:,cntr],Matx[:,-1])
    Matx_zscore=zscore(Matx[:,:72],axis=1)
    assign_vertex=np.argmax(Matx_zscore, axis=1)
    msort=np.sort(Matx_zscore,axis=1)
    mtokeep=msort[:,-1]-msort[:,-2]
    for avii in range(len(assign_vertex)):
        if Matx[avii,assign_vertex[avii]]<np.max(Matx[:,assign_vertex[avii]])/4.:
            assign_vertex[avii]=-100
            
    assign_vertex[np.logical_or(msort[:,-1]<3, mtokeep<=1)]=-100
    print "ctf finished!"
    
    
    if subject_inds[0] < 3:
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
# pick MEG channels
    picks = mne.pick_types(raw.info, meg=True, eeg=True, stim=False, eog=True, exclude='bads')

# Read epochs
    epochs = mne.Epochs(raw, events, event_ids, tmin, tmax, picks=picks, baseline=(None, 0), proj=True, reject=reject)
    epochs_concrete = epochs['Words']
    epochs_abstract = epochs['PWords']

    mne.epochs.equalize_epoch_counts([epochs_concrete,epochs_abstract])
    epochs=epochs_concrete
    events_no=epochs.copy().get_data().shape[0]
    

# First, we find the most active vertex in the left auditory cortex, which
#  we will later use as seed for the connectivity computation
#    snr = 3.0
#    lambda2 = 1.0 / snr ** 2
    evoked = epochs.average()
    stc_incoh = data_path + meg + 'firstMorphed_ico_typlex_ica_word_Source_Evoked_m500_700'
    stc_coh = mne.read_source_estimate(stc_incoh)
    #stc_coh.resample(100)

    for label1,lii in zip(labellist,label_index):
        print label1
        tmin, tmax = -0.5, 0.7
        if subject_inds[0] in np.array([3,4,13]):
            reject = dict(eeg=150e-6, grad=200e-12, mag=4e-12) #eog=150e-6,
        else:
            reject = dict(eeg=120e-6, grad=200e-12, mag=4e-12) #eog=150e-6,
        event_ids = {'cncrt_wrd': 1, 'abs_wrd': 2}
# frequency bands with variable number of cycles for wavelets
        frequencies = np.arange(8, 45, 1)
        n_cycles = frequencies / float(7)
        n_cycles[frequencies<=20] = 2
        n_cycles[frequencies<=20] += np.arange(0,1,0.077)
# n_cycles[np.logical_and((frequencies>10), (frequencies<=20))] = 3
        subject_from = subjects[subject_inds[0]]
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
        best_101thvert=sorted(bestvert_values,reverse=True)[1]
        best_100vert= best_verts[np.where(bestvert_values>best_101thvert)]
        print Matx_zscore[best_100vert,lii].mean()
        label_stc_data=np.zeros((stc_to.data.shape[0],1))
        label_stc_data[best_100vert]=1.
        label_stc=mne.SourceEstimate(label_stc_data, vertices=vertices_to,tmin=1e-3 * 1., tstep=1e-3 * 1., subject='fsaverage')
        if label1.name.endswith('-lh'):
            label_new=mne.stc_to_label(label_stc, smooth=False,subjects_dir=data_path,connected=False)[0]
            label_new.name=label1.name[:-3]+'_'+str(subject_no)+'_ctf-lh'
            thislabelname='lh'+label1.name[:-3]
        else:
            label_new=mne.stc_to_label(label_stc, smooth=False,subjects_dir=data_path,connected=False)[1]
            label_new.name=label1.name[:-3]+'_'+str(subject_no)+'_ctf-rh'
            thislabelname='rh'+label1.name[:-3]
        #label_new.save(label_new_path)
        print label_new.vertices.shape; print "label tc extracted"
        stcs=range(events_no)
        n_sources = stc_coh.data.shape[0]
        MI=np.zeros((n_sources,len(stcs),2))
        for ccc in range(len(stcs)):
            print ccc		
            data_in = data_path + meg + 'Morphed_ico_TypLex_epochs' + str(ccc) + '_icaclean_word_Source_Gamma_m500_700'
            stc_orig=mne.read_source_estimate(data_in)
            
            seed_ts = mne.extract_label_time_course(stc_orig, label_new, src, mode='mean_flip')
            for actc,tc in enumerate(range(150,251,100)): #time count
                print tc, actc
                curr_sig=stc_orig.data[:,tc+500:tc+500+200].copy()
                for vc in range(n_sources): #vertex count
                    MI[vc,ccc,actc]=mi(curr_sig[vc,:],seed_ts[0,tc+500:tc+500+200],bins=7)#[0,1]
                    
                
#            stc_orig.resample(100)
        MI_final=MI.mean(axis=1)
        
        tmin=0.15
        tstep=0.1
        
        MI_stc = mne.SourceEstimate(MI_final, vertices=stc_coh.vertices, tmin= tmin, tstep= tstep, subject='fsaverage')
        data_out = data_path + meg +  'Morphed_ico_typlex_ica_equalized_word_new_ctf_MI_Gamma_150_450_200ms_' + thislabelname                                        
        print data_out
        MI_stc.save(data_out)
        
        
