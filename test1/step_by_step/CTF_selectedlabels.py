# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 13:05:42 2016

@author: rf02
"""
print(__doc__)
print "hi! this script is for spectral functional connectivity. Very nice maps are fruits which worth waiting for I bet!"
import sys

#sys.path.append('/imaging/local/software/python_packages/nibabel/1.3.0')
#sys.path.append('/imaging/local/software/python_packages/pysurfer/v0.4')
#sys.path.append('/imaging/rf02/scikit-learn-0.15.0')
#
##sys.path.append('/imaging/local/software/python_packages/scikit-learn/v0.14.1')
#sys.path.append('/imaging/local/software/mne_python/la_v0.9')
sys.path.insert(1,'/imaging/rf02/TypLexMEG/mne_python0.9')
#sys.path.insert(1,'/imaging/local/software/mne_python/la_v0.9full')

sys.path.insert(1,'/imaging/rf02/scikit-learn-0.16.1')
sys.path.insert(1,'/home/rf02/.local/lib/python2.7/site-packages')


# for qsub
# (add ! here if needed) 
sys.path.append('/imaging/local/software/EPD/latest/x86_64/bin/python')
#sys.path.append('/imaging/rf02/TypLexMEG/mne_python_v8/')
#sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
###
import sklearn
from scipy.stats.mstats import zscore
import os
import numpy as np
import mne
from mne import io
from mne.io import Raw
import operator
from mne.minimum_norm import (read_inverse_operator, make_inverse_operator, write_inverse_operator, apply_inverse, apply_inverse_epochs,psf_ctf_new, point_spread_function,psf_ctf_18Dec15, psf_ctf_28Oct15)
from mne.minimum_norm.inverse import (prepare_inverse_operator)
import scipy.io
from mne import find_events
from mne.connectivity import seed_target_indices, spectral_connectivity
from mne import read_labels_from_annot
from mne.label import _n_colors as nc
from mne.label import _read_annot as mnera
from mne.label import _write_annot as mnewa
import copy
data_path = '/imaging/rf02/TypLexMEG/' # root directory for your MEG data
os.chdir(data_path)
event_path = '/group/erp/data/olaf.hauk/MEG/TypLexM/data/MF22_new/'
#label_path = '/group/erp/data/olaf.hauk/MEG/TypLexM/labels/createdlabels_SL/'
inv_path = '/imaging/rf02/TypLexMEG/'
subjects_dir = data_path 
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
labelss = mne.read_labels_from_annot(subject='fsaverage', parc='aparc', subjects_dir=data_path) #rez_aparc_24sep15
#labelss.remove(labelss[-3])
#labelss.remove(labelss[-2])
labelss.remove(labelss[-1])
for lbl in labelss:
    lbl.values.fill(1.0)
label_names=[label.name for label in labelss]
#label_index=np.array([label_names.index('ATL_latmed-lh'),label_names.index('anterior_inferiorfrontal-lh')])
#label_index=np.array([label_names.index('supramarginal-lh'),label_names.index('insula-lh')])
#label_index=np.array([label_names.index('parsorbitalis-lh'),label_names.index('parstriangularis-lh')])
label_index=np.array([label_names.index('postcentral-lh')])





labellist=[labelss[ii] for ii in label_index]
ii=-1
for cnt, meg in enumerate(list_all):
    
    print cnt, "preload variables"
    subjects_dir = data_path 
    fname_fwd = data_path + meg + 'ico5_forward_5-3L-EMEG-fwd.fif'
    fname_inv = inv_path + meg + 'InvOp_ico5newreg_fft_1_48_clean_ica_EMEG-inv.fif'
    vertices_avg = [np.arange(10242), np.arange(10242)]
    subject_no=subject_inds[cnt]

    forward = mne.read_forward_solution(fname_fwd,surf_ori=True, force_fixed=True)#, surf_ori=True)

    # read inverse operators
    inverse_operator_eegmeg = read_inverse_operator(fname_inv)
    
    snr = 3.0#ediv.mean()
    lambda2 = 1.0 / snr ** 2
    method = 'MNE'  # can be 'MNE' or 'sLORETA'
    mode = 'svd'
    n_svd_comp = 2
    lblname_cnt=0
    sub_src=inverse_operator_eegmeg['src']
    
    labels = copy.deepcopy(labellist)#mne.read_labels_from_annot(subject=subjects[subject_no], parc=subjects[subject_no]+'_morphed_split_aparc', subjects_dir=subjects_dir)#parc='rez_aparc_24sep15'
    #labels.remove(labels[-1])# unknown lh
    #labels.remove(labels[-1])# unknown rh
    labold=[label.name for label in labels]
    morphed_labels=[label.morph(subject_from='fsaverage', smooth=5, subject_to=subjects[subject_no], subjects_dir=data_path) for label in labels]
#        morphed_labels=list()
#        for jj,label in enumerate(labels):
#            print jj
#            morphed_labels.append(label.morph(subject_from='fsaverage', smooth=2, grade=vertices_avg, subject_to=subjects[subject_no], subjects_dir=data_path))
    n_svd_comp=1
    try:
        stc_ctf_mne=psf_ctf_28Oct15.cross_talk_function(inverse_operator_eegmeg, forward, labels=morphed_labels,method=method, lambda2=lambda2, signed=False, mode=mode, n_svd_comp=n_svd_comp, verbose=None)
        ii=ii+1
        print "psf/ctf finished!"
        subject_from = subjects[subject_no]
        subject_avg = 'fsaverage'
        morphmat1=mne.compute_morph_matrix(subject_from=subject_from, subject_to=subject_avg, vertices_from=stc_ctf_mne.vertices, vertices_to=vertices_avg, subjects_dir=data_path)
        stc_to=stc_ctf_mne.morph_precomputed(subject_avg,vertices_avg,morphmat1, subjects[subject_no])
        if ii==0:
            stc_all_data=stc_to.data[:,:,np.newaxis]
        else:
            stc_all_data=np.concatenate((stc_all_data,stc_to.data[:,:,np.newaxis]),axis=2)
    except:
        pass
    print stc_all_data.shape
    stc_sub=stc_all_data[:,0,:]-stc_all_data[:,1,:]
    stc_mean=stc_sub.mean(axis=1)
    
stc_summary=np.concatenate((stc_mean[:,np.newaxis],stc_all_data.mean(axis=2)),axis=1)
summary_stc = mne.SourceEstimate(stc_summary, vertices=vertices_avg, tmin=1e-3*2, tstep=1e-3*2, subject=subject_avg, verbose=None)
out_file=data_path + 'postcentral_CTF'
summary_stc.save(out_file)
import matplotlib.pyplot as plt
n,b,p= plt.hist(stc_all_data.mean(axis=2)[:,2],bins=100)
plt.xlabel('CTF value')
plt.ylabel('Number of Vertices')
plt.title('Supramarginal AND LateralOccipital')
    