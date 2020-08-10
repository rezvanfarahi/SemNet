# -*- coding: utf-8 -*-
"""
Created on Sat Dec 19 17:43:49 2015

@author: rf02
"""

"""
==========================================================
Compute point-spread functions (PSFs) for MNE/dSPM/sLORETA
==========================================================

PSFs are computed for four labels in the MNE sample data set
for linear inverse operators (MNE, dSPM, sLORETA).
PSFs describe the spread of activation from one label
across the cortical surface.
"""

# Authors: Olaf Hauk <olaf.hauk@mrc-cbu.cam.ac.uk>
#          Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
#
# License: BSD (3-clause)



# Author: Martin Luessi <mluessi@nmr.mgh.harvard.edu>
#
# License: BSD (3-clause)

print(__doc__)
print "hi! this script is for spectral functional connectivity. Very nice maps are fruits which worth waiting for I bet!"
import sys

#sys.path.append('/imaging/local/software/python_packages/nibabel/1.3.0')
#sys.path.append('/imaging/local/software/python_packages/pysurfer/v0.4')
#sys.path.append('/imaging/rf02/scikit-learn-0.15.0')
#
##sys.path.append('/imaging/local/software/python_packages/scikit-learn/v0.14.1')
#sys.path.append('/imaging/local/software/mne_python/latest_v0.9')
sys.path.insert(1,'/imaging/rf02/TypLexMEG/mne_python0.9')
#sys.path.insert(1,'/imaging/local/software/mne_python/latest_v0.9full')

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
from mne.minimum_norm import (read_inverse_operator, make_inverse_operator, write_inverse_operator, apply_inverse, apply_inverse_epochs,psf_ctf_new, point_spread_function,psf_ctf_18Dec15)
from mne.minimum_norm.inverse import (prepare_inverse_operator)
import scipy.io
from mne import find_events
from mne.connectivity import seed_target_indices, spectral_connectivity
from mne import read_labels_from_annot
from mne.label import _n_colors as nc
from mne.label import _read_annot as mnera
from mne.label import _write_annot as mnewa
from mne.label import write_labels_to_annot 
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

ll = []
for ss in subject_inds:
 ll.append(list_all[ss])


print "ll:"
print ll
ii=-1
labelss = mne.read_labels_from_annot(subject='fsaverage', parc='aparc.a2009s', subjects_dir=subjects_dir)#parc='rez_aparc_24sep15'
#labelss.remove(labelss[-3])
#labelss.remove(labelss[-2])
labelss.remove(labelss[-1])
labelss.remove(labelss[-1])

print len(labelss)
llk=labelss[0]+labelss[-2]
lrk=labelss[1]+labelss[-1]
llk.name=labelss[0].name
lrk.name=labelss[1].name
#labelss.remove(labelss[-1])
#labelss.remove(labelss[-1])
Matxl=np.zeros((20484,17))
Matxr=np.zeros((20484,17))
Matx=np.zeros((20484,len(labelss)+1,17))
rat=np.zeros((len(labelss),len(labelss),17))
rat_blnc=np.zeros((len(labelss),17))
stc_all=range(17)
overlap_sub_winner=np.zeros((len(labelss),17))
assign_vertex=np.zeros((20484,17))
exp90=np.zeros((len(labelss),17))
label_colors = [label.color for label in labelss]
# read label(s)
for lbl in labelss:
    lbl.values.fill(1.0)
for cnt, meg in enumerate(ll):
    ## load variables
    ii=ii+1
    print cnt
    subjects_dir = data_path 
    fname_fwd = data_path + meg + 'ico5_forward_5-3L-EMEG-fwd.fif'
    fname_inv = inv_path + meg + 'InvOp_ico5newreg_fft_1_48_clean_ica_EMEG-inv.fif'
    vertices_to = [np.arange(10242), np.arange(10242)]
    subject_no=subject_inds[cnt]
    
    forward = mne.read_forward_solution(fname_fwd,surf_ori=True, force_fixed=True)#, surf_ori=True)

    # read inverse operators
    inverse_operator_eegmeg = read_inverse_operator(fname_inv)
    sub_src=inverse_operator_eegmeg['src']    
    ## morph lables
    labels = [labelss[jj].morph(subject_from='fsaverage', smooth=5, subject_to=subjects[subject_no], subjects_dir=data_path) for jj in range(len(labelss))]
    print "label morph fin"
    labels_forcolor = copy.deepcopy(labels)
    subject_from = subjects[subject_no]
    subject_to = 'fsaverage'
    ## arange for new colorings etc.
    annotl=np.zeros((sub_src[0]['np'],1))
    annotr=np.zeros((sub_src[1]['np'],1))
    CTAB_L=np.zeros((len(labels)/2,5))
    CTAB_R=np.zeros((len(labels)/2,5))
    colors=nc(len(labels)/2,bytes_=True, cmap='hsv')
    colors=colors[np.random.permutation(len(labels)/2)]
    annot_id_coding = np.array((1, 2 ** 8, 2 ** 16))
    annot_id=np.sum(colors[:,:3] * annot_id_coding, axis=1)
    for lcntn in range(len(labels)):
        print lcntn
        if labels[lcntn].name.endswith('lh'):
            annotl[labels[lcntn].vertices]=annot_id[lcntn/2]
        else:
            annotr[labels[lcntn].vertices]=annot_id[int(lcntn/2)]
    CTAB_L[:,:4]=colors; CTAB_L[:,4]=annot_id
    CTAB_R[:,:4]=colors; CTAB_R[:,4]=annot_id
    label_names=[labels[lh].name[:-3] for lh in range(len(labels))]
    NAMES_L=label_names[::2]
    NAMES_R=label_names[1::2]
    ## write morphed labels
    import pickle
    fnamem=data_path + subject_from+'/label/morphed_aparc_2009_09Jul16'
    with open(fnamem,'wb') as f:
        pickle.dump(labels,f)
    
    fnamel=data_path + subject_from+'/label/lh.morphed_aparc_2009_21Dec15.annot'
    mnewa(fnamel, np.squeeze(annotl), CTAB_L, NAMES_L)

    fnamer=data_path + subject_from+'/label/rh.morphed_aparc_2009_21Dec15.annot'
    mnewa(fnamer, np.squeeze(annotr), CTAB_R, NAMES_R)
    
    annotlorig,ctablorig,nameslorig=mnera('/imaging/rf02/TypLexMEG/fsaverage/label/lh.aparc.a2009s.annot')
    annotrorig,ctabrorig,namesrorig=mnera('/imaging/rf02/TypLexMEG/fsaverage/label/rh.aparc.a2009s.annot')
    ## rewrite original labels with the same colors
    fnameorigl=data_path +'fsaverage/label/lh.aparc_2009_21Dec15.annot'
    mnewa(fnameorigl, np.squeeze(annotlorig), CTAB_L, NAMES_L)

    fnameorigr=data_path + 'fsaverage/label/rh.aparc_2009_21Dec15.annot'
    mnewa(fnameorigr, np.squeeze(annotrorig), CTAB_R, NAMES_R)
    
    ## splitting 
    snr = 3.0#ediv.mean()
    lambda2 = 1.0 / snr ** 2
    method = 'MNE'  # can be 'MNE' or 'sLORETA'
    mode = 'svd'
    n_svd_comp = 1
    lblname_cnt=0
    stc_ctf_mne, explained90=psf_ctf_18Dec15.cross_talk_function(inverse_operator_eegmeg, forward, labels,method=method, lambda2=lambda2, signed=False, mode=mode, n_svd_comp=1, verbose=None)
    exp90[:,cnt]=np.squeeze(explained90)#.astype(int)
EXP90p=np.zeros((exp90.shape[0],1))
EXP90n=np.zeros((exp90.shape[0],1))
for cnt in range(exp90.shape[0]):
    EXP90p[cnt]=scipy.stats.mstats.mode(exp90[cnt,:].astype(int))[0][0]
    EXP90n[cnt]=scipy.stats.mstats.mode(exp90[cnt,:].astype(int))[1][0]
EXP90p=np.squeeze(EXP90p)
for c90 in range(int(len(labels)/2)):
    pmin=np.min([EXP90p[c90*2],EXP90p[c90*2+1]])
    EXP90p[c90*2]=pmin; EXP90p[c90*2+1]=pmin
EXP90=EXP90p+1
path_c=data_path+'explained90all_2009.mat'
scipy.io.savemat(path_c,{'EXP90':EXP90})
#            
#    llnew=[labels[ii] for ii in range(len(labels)) if explained90[ii][0]==1.]
#    lln1=[labels[ii].split(explained90[ii][0], subjects_dir=data_path)[:] for ii in range(len(labels)) if explained90[ii][0]>1]
#    llnn=list()
#    for ii in range(len(lln1)):
#        llnn+=lln1[ii]
#    print "psf/ctf finished!"
#    labels_backup=copy.deepcopy(labels)
#    labels=llnn+llnew
#    stc_ctf_mne, explained90=psf_ctf_18Dec15.cross_talk_function(inverse_operator_eegmeg, forward, labels,method=method, lambda2=lambda2, signed=False, mode=mode, n_svd_comp=1, verbose=None)
    

#    rs=nc(len(labels),bytes_=False, cmap='hsv')
#    rs=rs[np.random.permutation(len(labels))]
#    import copy
#    labels_removeov=copy.deepcopy((labels))
#    ovt=np.array([])
#    for cnt1 in range(len(labels)):
#        print cnt1
#        for cnt2 in range(len(labels)):
#            if cnt1<cnt2:
#                label1ver=labels_removeov[cnt1].vertices
#                label2ver=labels_removeov[cnt2].vertices
#                vnoverlap2=np.in1d(label2ver,label1ver,invert=True)#voverlap=np.intersect1d(label1ver,label2ver)
#                label2ver_new=label2ver[vnoverlap2]#np.delete(label2ver,overlap_index)
#                jjh=labels_removeov[cnt2]
#                jjh.vertices=label2ver_new
#                vnoverlap1=np.in1d(label1ver,label2ver,invert=True)
#                voverlap=np.in1d(label1ver,label2ver)
#                label1ver_new=label1ver[vnoverlap1]#np.delete(label2ver,overlap_index)
#                ovver_new=label1ver[voverlap]
#                jji=labels_removeov[cnt1]
#                jji.vertices=label1ver_new
#                ovt=np.concatenate((ovver_new,ovt),axis=0)
#        labels_removeov[cnt1].color=labels[cnt1].color
#    labels_semifinal=[labels_removeov[ii] for ii in range(len(labels))]
#    ov_label=copy.deepcopy(labels_semifinal[0])
#    ov_label.vertices=ovt
#    ov_label.name='overlaps'
#    ov_label.color=rs[0]
#    labels_semifinal.append(ov_label)
#    mne.write_labels_to_annot(labels_semifinal, subject=subject_from, parc=subject_from+'_first_morphed_parc', subjects_dir=data_path,colormap='hsv',hemi='both',overwrite=True)
#    