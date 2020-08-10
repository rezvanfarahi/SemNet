# -*- coding: utf-8 -*-
"""
Created on Wed Dec 23 14:34:39 2015

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
from mne.minimum_norm import (read_inverse_operator, make_inverse_operator, write_inverse_operator, apply_inverse, apply_inverse_epochs,psf_ctf_new, point_spread_function,psf_ctf_18Dec15)
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

ll = []
for ss in subject_inds:
 ll.append(list_all[ss])


print "ll:"
print ll
ii=-1

labels_fc = mne.read_labels_from_annot(subject='fsaverage', parc='fsaverage_morphed_split_aparc', subjects_dir=subjects_dir)#parc='rez_aparc_24sep15'
labels_fc.remove(labels_fc[-1])# unknown lh
labels_fc.remove(labels_fc[-1])# unknown rh
#
split_lab_verts=[label.vertices.shape[0] for label in labels_fc]
min_act_lbl=np.min(split_lab_verts)
Matx=np.zeros((20484,len(labels_fc)+1,17))

stc_all=range(17)
assign_vertex=np.zeros((20484,17))

label_colors = [label.color for label in labels_fc]
# read label(s)
for lbl in labels_fc:
    lbl.values.fill(1.0)
for cnt, meg in enumerate(ll):
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
    labels = mne.read_labels_from_annot(subject=subjects[subject_no], parc=subjects[subject_no]+'_morphed_split_aparc', subjects_dir=subjects_dir)#parc='rez_aparc_24sep15'
    labels.remove(labels[-1])# unknown lh
    labels.remove(labels[-1])# unknown rh
    label_names=[labels[lcn3].name for lcn3 in range(len(labels))]
    
    #labels = [labelss[jj].morph(subject_from='fsaverage', smooth=5, subject_to=subjects[subject_no], subjects_dir=data_path) for jj in range(len(labelss))]
    labels_forcolor = copy.deepcopy(labels)

    snr = 3.0#ediv.mean()
    lambda2 = 1.0 / snr ** 2
    method = 'MNE'  # can be 'MNE' or 'sLORETA'
    mode = 'svd'
    n_svd_comp = 1
    lblname_cnt=0
    sub_src=inverse_operator_eegmeg['src']
    
    ## make and write divided annot
    
    subject_from = subjects[subject_no]
    subject_to = 'fsaverage'
    colors=nc(len(labels),bytes_=False, cmap='gist_ncar')#Note! left right different colors, not important here.
    colors=colors[np.random.permutation(int(len(labels)))]
    
    stc_ctf_mne, ep80=psf_ctf_18Dec15.cross_talk_function(inverse_operator_eegmeg, forward, labels,method=method, lambda2=lambda2, signed=False, mode=mode, n_svd_comp=1, verbose=None)
    lv=[len(labels[vc].vertices) for vc in range(len(labels))]
    #mne.write_labels_to_annot(labels, subject=subject_from, parc=subject_from+'_preparc_smooth_split', subjects_dir=data_path,colormap='gist_ncar',hemi='both',overwrite=True)
    morphmat1=mne.compute_morph_matrix(subject_from, subject_to, stc_ctf_mne.vertices, vertices_to, subjects_dir=data_path)
    stc_to=stc_ctf_mne.morph_precomputed(subject_to,vertices_to,morphmat1, subject_from)
    if cnt==0:
        stc_all_data=stc_to.data[:,:,np.newaxis]
    else:
        stc_all_data=np.concatenate((stc_all_data,stc_to.data[:,:,np.newaxis]),axis=2)

vertices_avg= [np.arange(10242), np.arange(10242)]
subject_avg='fsaverage'
srcin='/imaging/rf02/TypLexMEG/fsaverage/bem/fsaverage_dist-ico-5-src.fif'
#src_avg2=src_avg.copy()
src_avg = mne.read_source_spaces(srcin)
mne.add_source_space_distances(src_avg)
stc_mean_data=np.mean(stc_all_data,axis=2)
Matx=stc_mean_data[:,:-1]# Nv x Nlabels+1
path_mtx=data_path+'Matx_first.mat'
scipy.io.savemat(path_mtx,{'stc_mean_data':stc_mean_data})
Matx_zscore=np.zeros(Matx.shape)
#    for cntr in range(Matx.shape[1]-1):
#        Matx_normalised[:,cntr]=np.divide(Matx[:,cntr],Matx[:,-1])
Matx_zscore=zscore(Matx,axis=1)
assign_vertex=np.argmax(Matx_zscore, axis=1)
msort=np.sort(Matx_zscore,axis=1)
mtokeep=msort[:,-1]-msort[:,-2]
for avii in range(len(assign_vertex)):
    if Matx[avii,assign_vertex[avii]]<np.max(Matx[:,assign_vertex[avii]])/4.:
        assign_vertex[avii]=-100
        
assign_vertex[np.logical_or(msort[:,-1]<3, mtokeep<=1)]=-100
uv=np.unique(assign_vertex)
#num_oldv=np.zeros((len(uv)-1,1))
#for olcnt in range(len(uv)-1):
#    num_oldv[olcnt]=np.where(assign_vertex==uv[olcnt+1])[0].shape[0]
##make merged labels
merge_candidates=np.where(np.logical_and(np.logical_and(msort[:,-1]>=3,msort[:,-2]>=3), mtokeep<=1))[0]
margsort=np.argsort(Matx_zscore,axis=1)
merge_labels1=margsort[merge_candidates,-2:]
merge_labels=np.sort(merge_labels1,axis=1)
new_labels=[tuple(row) for row in merge_labels]
newlabels=np.unique(new_labels)
labels_merged=list()
lblname_cnt=-1
for labcnt in range(newlabels.shape[0]):
    print labcnt
    thislabel_verts=np.where(np.logical_and(merge_labels[:,0]==newlabels[labcnt,0],merge_labels[:,1]==newlabels[labcnt,1]))[0]
    thislabel=np.asarray([merge_candidates[vrt] for vrt in thislabel_verts])
    print thislabel.shape
    lblname_cnt=lblname_cnt+1
    thislabel_data=np.zeros((assign_vertex.shape[0],1))
    thislabel_data[thislabel]=1.0
    thislabel_stc = mne.SourceEstimate(thislabel_data, vertices=vertices_avg, tmin=1e-3*2, tstep=1e-3*2, subject=subject_avg, verbose=None)
    label_append=[x for x in mne.stc_to_label(thislabel_stc,src=src_avg, smooth=True,connected=False, subjects_dir=data_path) if x is not None]
    if thislabel.shape[0]>=min_act_lbl:
        labels_merged.append(label_append)
        labels_merged[lblname_cnt][0].name=labels[newlabels[labcnt,1]].name[:-3]+'_'+labels[newlabels[labcnt,0]].name
unver=np.unique(assign_vertex)

labvernum=np.zeros(unver.shape[0])
for uc in range(len(unver)):
    labvernum[uc]=np.where(assign_vertex==unver[uc])[0].shape[0]
    if labvernum[uc]<20.0:
        assign_vertex[np.where(assign_vertex==unver[uc])[0]]=-100
unver=np.unique(assign_vertex)
print unver
labvernum=np.zeros(unver.shape[0])
for uc in range(len(unver)):
    labvernum[uc]=np.where(assign_vertex==unver[uc])[0].shape[0]        
print sorted(labvernum)
print labvernum.shape
for lcnt in range(len(unver)-1):
    thislabel=np.where(assign_vertex==unver[lcnt+1])[0]
    thislabel_data=np.zeros((assign_vertex.shape[0],1))
    thislabel_data[thislabel]=1.0
    thislabel_stc = mne.SourceEstimate(thislabel_data, vertices=vertices_avg, tmin=1e-3*2, tstep=1e-3*2, subject=subject_avg, verbose=None)
    label_append=[x for x in mne.stc_to_label(thislabel_stc,src=src_avg, smooth=True,connected=False, subjects_dir=data_path) if x is not None]
    labels_merged.append(label_append)
    labels_merged[lblname_cnt+lcnt+1][0].name=labels[unver[lcnt+1]].name
labelss=[label[0] for label in labels_merged]
labelss_names=[label.name for label in labelss]
balanced_labels01=list()
for lcnt01 in range(len(labelss)):
    if labelss_names[lcnt01].endswith('lh'):
        for lcnt02 in range(len(labelss)):
            if labelss_names[lcnt02].endswith('rh') and labelss_names[lcnt02][:-3]==labelss_names[lcnt01][:-3]:
                balanced_labels01.append(labelss[lcnt01])
                balanced_labels01.append(labelss[lcnt02])

import copy
#rs=nc(len(labels_final),bytes_=False, cmap='gist_ncar')
#rs=rs[np.random.permutation(len(labels_final))]
annotl=np.zeros((src_avg[0]['np'],1))
annotr=np.zeros((src_avg[1]['np'],1))
CTAB_L=np.zeros((len(balanced_labels01)/2,5))
CTAB_R=np.zeros((len(balanced_labels01)/2,5))
colors=nc(int(len(balanced_labels01)/2),bytes_=True, cmap='gist_ncar')
colors=colors[np.random.permutation(int(len(balanced_labels01)/2))]
annot_id_coding = np.array((1, 2 ** 8, 2 ** 16))
annot_id=np.sum(colors[:,:3] * annot_id_coding, axis=1)
for lcntn in range(len(balanced_labels01)):
    print lcntn
    if balanced_labels01[lcntn].name.endswith('lh'):
        annotl[balanced_labels01[lcntn].vertices]=annot_id[lcntn/2]
    else:
        annotr[balanced_labels01[lcntn].vertices]=annot_id[int(lcntn/2)]
CTAB_L[:,:4]=colors; CTAB_L[:,4]=annot_id
CTAB_R[:,:4]=colors; CTAB_R[:,4]=annot_id
label_names=[balanced_labels01[lh].name[:-3] for lh in range(len(balanced_labels01))]
NAMES_L=label_names[::2]
NAMES_R=label_names[1::2]
## write morphed labels
fnamel=data_path + subject_avg+'/label/lh.splitlabels_merged_labels_prefinal.annot'
mnewa(fnamel, np.squeeze(annotl), CTAB_L, NAMES_L)

fnamer=data_path + subject_avg+'/label/rh.splitlabels_merged_labels_prefinal.annot'
mnewa(fnamer, np.squeeze(annotr), CTAB_R, NAMES_R)

for cnt, meg in enumerate(ll):
    ii=ii+1
    print cnt
    subjects_dir = data_path 
    fname_fwd = data_path + meg + 'ico5_forward_5-3L-EMEG-fwd.fif'
    fname_inv = inv_path + meg + 'InvOp_ico5newreg_fft_1_48_clean_ica_EMEG-inv.fif'
    vertices_avg = [np.arange(10242), np.arange(10242)]
    subject_no=subject_inds[cnt]

    forward = mne.read_forward_solution(fname_fwd,surf_ori=True, force_fixed=True)#, surf_ori=True)

    # read inverse operators
    inverse_operator_eegmeg = read_inverse_operator(fname_inv)
    
    labels = [balanced_labels01[jj].morph(subject_from='fsaverage', smooth=5, subject_to=subjects[subject_no], subjects_dir=data_path) for jj in range(len(balanced_labels01))]
    
    snr = 3.0#ediv.mean()
    lambda2 = 1.0 / snr ** 2
    method = 'MNE'  # can be 'MNE' or 'sLORETA'
    mode = 'svd'
    n_svd_comp = 1
    lblname_cnt=0
    sub_src=inverse_operator_eegmeg['src']
    
    #final ctf
    stc_ctf_mne,expl80=psf_ctf_18Dec15.cross_talk_function(inverse_operator_eegmeg, forward, labels,method=method, lambda2=lambda2, signed=False, mode=mode, n_svd_comp=1, verbose=None)
    print "psf/ctf finished!"
    morphmat1=mne.compute_morph_matrix(subject_from=subjects[subject_no], subject_to='fsaverage', vertices_from=stc_ctf_mne.vertices, vertices_to=vertices_avg, subjects_dir=data_path)
    stc_to=stc_ctf_mne.morph_precomputed(subject_avg,vertices_avg,morphmat1, subjects[subject_no])
    if cnt==0:
        stc_all_data2=stc_to.data[:,:,np.newaxis]
    else:
        stc_all_data2=np.concatenate((stc_all_data2,stc_to.data[:,:,np.newaxis]),axis=2)
stc_mean_data2=np.mean(stc_all_data2,2)
Matx=stc_mean_data2[:,:-1]# Nv x Nlabels+1
path_mtx=data_path+'Matx_second.mat'
scipy.io.savemat(path_mtx,{'stc_mean_data':stc_mean_data})
Matx_zscore=np.zeros(Matx.shape)
Matx_zscore=zscore(Matx,axis=1)
assign_vertex=np.argmax(Matx_zscore, axis=1)
msort=np.sort(Matx_zscore,axis=1)
mtokeep=msort[:,-1]-msort[:,-2]
for avii in range(len(assign_vertex)):
    if Matx[avii,assign_vertex[avii]]<np.max(Matx[:,assign_vertex[avii]])/4.:
        assign_vertex[avii]=-100
assign_vertex[np.logical_or(msort[:,-1]<3, mtokeep<=1)]=-100
unver=np.unique(assign_vertex)
print unver
labvernum=np.zeros(unver.shape[0])
for uc in range(len(unver)):
    labvernum[uc]=np.where(assign_vertex==unver[uc])[0].shape[0] 
print sorted(labvernum)
print labvernum.shape


colovbal=list()
colovfc=list()
for cc1 in range(len(labels_fc)):
    for cc2 in range(len(balanced_labels01)):
        if labels_fc[cc1].name==balanced_labels01[cc2].name:
            colovbal.append((cc2))
            colovfc.append((cc1))
            labels_fc[cc1].color=balanced_labels01[cc2].color
        
fcnov=np.arange(len(labels_fc))[np.in1d(np.arange(len(labels_fc)),np.asarray(colovfc),invert=True)]
colors2=nc(int((len(balanced_labels01)+len(labels_fc))/2),bytes_=False, cmap='gist_ncar')[int(len(balanced_labels01)/2)+1:]
colors2=colors2[np.random.permutation(len(colors2))]
for cc2,cc1 in enumerate(fcnov):
    labels_fc[cc1].color=colors2[int(cc2/2),:]
            
#mne.write_labels_to_annot(labels_fc, subject=subject_avg, parc=subject_avg+'_preparc_split_synccolor', subjects_dir=data_path,colormap='gist_ncar',hemi='both',overwrite=True)
labelsf=list()
for flc in range(len(unver)-1):
    print flc
    thislabel=np.where(assign_vertex==unver[flc+1])[0]
    thislabel_data=np.zeros((assign_vertex.shape[0],1))
    thislabel_data[thislabel]=1.0
    thislabel_stc = mne.SourceEstimate(thislabel_data, vertices=vertices_avg, tmin=1e-3*2, tstep=1e-3*2, subject=subject_avg, verbose=None)
    label_append=[x for x in mne.stc_to_label(thislabel_stc,src=src_avg, smooth=True,connected=False, subjects_dir=data_path) if x is not None]
    labelsf.append(label_append)
    labelsf[flc][0].name=balanced_labels01[unver[flc+1]].name
    labelsf[flc][0].color=balanced_labels01[unver[flc+1]].color
labelsff=[label[0] for label in labelsf]
lnames=[labelsff[lcnt01].name for lcnt01 in range(len(labelsff))]
balanced_labels=list()
for lcnt01 in range(len(labelsff)):
    if lnames[lcnt01].endswith('lh'):
        for lcnt02 in range(len(labelsff)):
            if lnames[lcnt02].endswith('rh') and lnames[lcnt02][:-3]==lnames[lcnt01][:-3]:
                balanced_labels.append(labelsff[lcnt01])
                balanced_labels.append(labelsff[lcnt02])
                
balanced_labels02 = mne.read_labels_from_annot(subject=subject_avg, parc='splitlabels_merged_labels_prefinal', subjects_dir=subjects_dir)
for iii in range(len(balanced_labels02)):
    for jjj in range(len(balanced_labels)):
        if balanced_labels[jjj].name==balanced_labels02[iii].name:
            balanced_labels[jjj].color=balanced_labels02[iii].color
        
mne.write_labels_to_annot(balanced_labels, subject=subject_avg, parc=subject_avg+'_finalpostparc_smooth_split_merge', subjects_dir=data_path,colormap='gist_ncar',hemi='both',overwrite=True)

#path_c=data_path+'colormatrix.mat'
#scipy.io.savemat(path_c,{'rs':rs})
#lfname=[labelsff[ii].name for ii in range(len(labelsff))]
#lfname_path=data_path+subject_from+'/labels_modified_names_split_merge.mat'
#scipy.io.savemat(lfname_path,{'lfname':lfname})
