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
from mne.minimum_norm import (read_inverse_operator, make_inverse_operator, write_inverse_operator, apply_inverse, apply_inverse_epochs, point_spread_function,psf_ctf_18Dec15, psf_ctf_28Oct15)
from mne.minimum_norm.inverse import (prepare_inverse_operator)
import scipy.io
from mne import find_events
from mne.connectivity import seed_target_indices, spectral_connectivity
from mne import read_labels_from_annot
from mne.label import _n_colors as nc
from mne.label import _read_annot as mnera
from mne.label import _write_annot as mnewa
import copy
import pickle
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

def make_parcellation_mine(src,labels,colors,subject,fname):
    annotl=np.zeros((src[0]['np'],1))
    annotr=np.zeros((src[1]['np'],1))
    CTAB_L=np.zeros((len(labels)/2,5))
    CTAB_R=np.zeros((len(labels)/2,5))
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
    fnamel=data_path + subject+'/label/lh.'+subject+'_'+fname+'.annot'
    mnewa(fnamel, np.squeeze(annotl), CTAB_L, NAMES_L)
    
    fnamer=data_path + subject+'/label/rh.'+subject+'_'+fname+'.annot'
    mnewa(fnamer, np.squeeze(annotr), CTAB_R, NAMES_R)

ll = []
for ss in subject_inds:
 ll.append(list_all[ss])


print "ll:"
print ll
ii=-1

labels_fc = mne.read_labels_from_annot(subject='fsaverage', parc='fsaverage_morphed_split_aparc_2009', subjects_dir=subjects_dir)#parc='rez_aparc_24sep15'
#labels_test = mne.read_labels_from_annot(subject='fsaverage', parc='fsaverage_onlyvertices_actualRG_82', subjects_dir=subjects_dir)
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
    fnamemf=data_path + subjects[subject_no]+'/label/'+subjects[subject_no]+'_morphed_split_aparc_2009'
    with open(fnamemf, 'rb') as f:
        labels=pickle.load(f)
    #    labels = mne.read_labels_from_annot(subject=subjects[subject_no], parc=subjects[subject_no]+'_morphed_split_aparc', subjects_dir=subjects_dir)#parc='rez_aparc_24sep15'
    #    labels.remove(labels[-1])# unknown lh
    #    labels.remove(labels[-1])# unknown rh
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
#    colors=nc(len(labels),bytes_=False, cmap='gist_ncar')#Note! left right different colors, not important here.
#    colors=colors[np.random.permutation(int(len(labels)))]
    
    
    try: 
        stc_ctf_mne=psf_ctf_28Oct15.cross_talk_function(inverse_operator_eegmeg, forward, labels,method=method, lambda2=lambda2, signed=False, mode=mode, n_svd_comp=1, verbose=None)
        
        #mne.write_labels_to_annot(labels, subject=subject_from, parc=subject_from+'_preparc_smooth_split', subjects_dir=data_path,colormap='gist_ncar',hemi='both',overwrite=True)
        morphmat1=mne.compute_morph_matrix(subject_from, subject_to, stc_ctf_mne.vertices, vertices_to, subjects_dir=data_path)
        stc_to=stc_ctf_mne.morph_precomputed(subject_to,vertices_to,morphmat1, subject_from)
        if cnt==0:
            stc_all_data=stc_to.data[:,:,np.newaxis]
        else:
            stc_all_data=np.concatenate((stc_all_data,stc_to.data[:,:,np.newaxis]),axis=2)
    except:
        pass

print "stage 1 finished- ctf of split labels"
vertices_avg= [np.arange(10242), np.arange(10242)]
subject_avg='fsaverage'
srcin='/imaging/rf02/TypLexMEG/fsaverage/bem/fsaverage_dist-ico-5-src.fif'
#src_avg2=src_avg.copy()
src_avg = mne.read_source_spaces(srcin)
#mne.add_source_space_distances(src_avg)
stc_mean_data=np.mean(stc_all_data,axis=2)
Matx=stc_mean_data[:,:-1]# Nv x Nlabels+1
path_mtx=data_path+'rev1_Matx_first2009.mat'
scipy.io.savemat(path_mtx,{'stc_mean_data':stc_mean_data})
#ltc=scipy.io.loadmat(path_mtx,mat_dtype=True); stc_mean_data=np.squeeze(ltc['stc_mean_data'])
#Matx=stc_mean_data
Matx_zscore=np.zeros(Matx.shape)
#    for cntr in range(Matx.shape[1]-1):
#        Matx_normalised[:,cntr]=np.divide(Matx[:,cntr],Matx[:,-1])
Matx_zscore=zscore(Matx,axis=1)
assign_vertex=np.argmax(Matx_zscore, axis=1)
msort=np.sort(Matx_zscore,axis=1)
mtokeep=msort[:,-1]-msort[:,-2]
for avii in range(len(assign_vertex)):
    if Matx[avii,assign_vertex[avii]]<np.max(Matx[:,assign_vertex[avii]])/2.:
        assign_vertex[avii]=-100
        
assign_vertex[np.logical_or(msort[:,-1]<3, mtokeep<=1.0)]=-100
uv=np.unique(assign_vertex)
#num_oldv=np.zeros((len(uv)-1,1))
#for olcnt in range(len(uv)-1):
#    num_oldv[olcnt]=np.where(assign_vertex==uv[olcnt+1])[0].shape[0]

##make merged labels
merge_candidates=np.where(np.logical_and(np.logical_and(msort[:,-1]>=3,msort[:,-2]>=3), mtokeep<=1.0))[0]
margsort=np.argsort(Matx_zscore,axis=1)
merge_labels1=margsort[merge_candidates,-2:]

merge_labels=np.sort(merge_labels1,axis=1)
new_labels=[tuple(row) for row in merge_labels]
#newlabels=np.unique(new_labels)
newlabels=[]
for new_label in new_labels:
    if new_label not in newlabels:
        newlabels.append(new_label)
newlabels=np.asarray(newlabels)
labels_merged=list()
lblname_cnt=-1
nver_merge=list()
for labcnt in range(newlabels.shape[0]):
    print labcnt
    thislabel_verts=np.where(np.logical_and(merge_labels[:,0]==newlabels[labcnt,0],merge_labels[:,1]==newlabels[labcnt,1]))[0]
    thislabel=np.asarray([merge_candidates[vrt] for vrt in thislabel_verts])    
    thislabel_data=np.zeros((assign_vertex.shape[0],1))
    thislabel_data[thislabel]=1.0
    thislabel_stc = mne.SourceEstimate(thislabel_data, vertices=vertices_avg, tmin=1e-3*2, tstep=1e-3*2, subject=subject_avg, verbose=None)
    label_append=[x for x in mne.stc_to_label(thislabel_stc,src=src_avg, smooth=True,connected=False, subjects_dir=data_path) if x is not None]
    if len(label_append)>1:
        labh=np.array([len(label_append[0].vertices),len(label_append[1].vertices)])
        lain=np.argmax(labh)
    else:
        lain=0
    if label_append[lain].vertices.shape[0]>=min_act_lbl and labels[newlabels[labcnt,1]].hemi==labels[newlabels[labcnt,0]].hemi and label_append[lain].hemi==labels[newlabels[labcnt,0]].hemi:
        nver_merge.append(thislabel)
        lblname_cnt=lblname_cnt+1
        labels_merged.append(label_append[0])
        labels_merged[lblname_cnt].name=labels[newlabels[labcnt,1]].name[:-3]+'_'+labels[newlabels[labcnt,0]].name
unver=np.unique(assign_vertex)
#remove merged label vertices from the original label vertices
import copy
labels_split_merge=labels_fc+labels_merged
split_names=[label.name for label in labels_fc]
import copy
labels_split_merge=labels_fc+labels_merged
split_names=[label.name for label in labels_fc]
#path_mtx='/imaging/rf02/parcellation/MNIlh.mat'
#MNIlh_pre=scipy.io.loadmat(path_mtx,mat_dtype=True); MNIlh=np.squeeze(MNIlh_pre['MNIlh'])
#path_mtx='/imaging/rf02/parcellation/MNIrh.mat'
#MNIrh_pre=scipy.io.loadmat(path_mtx,mat_dtype=True); MNIrh=np.squeeze(MNIrh_pre['MNIrh'])
#
#nv_merge_lh=np.asarray([len(label.vertices) for label in labels_merged if label.name.endswith('lh')]).sum()
#nv_merge_rh=np.asarray([len(label.vertices) for label in labels_merged if label.name.endswith('rh')]).sum()
#
##labels_lh_ind=[ii1 for ii1 in range(len(unver)) if unver[ii1] in labels_lh_ind]
#
#
#labels_lh_ind=[llcnt for llcnt in range(len(labels_merged)) if labels_merged[llcnt].name.endswith('lh')]
#labels_rh_ind=[rrcnt for rrcnt in range(len(labels_merged)) if labels_merged[rrcnt].name.endswith('rh')]
##nver_merge=np.asarray(nver_merge)
#
#if nv_merge_lh>=nv_merge_rh:#left wins
#    unver=[nver_merge[lhi] for lhi in labels_lh_ind]
#    whichhem=1
#    labels_bh_ind=copy.deepcopy(labels_lh_ind)
#else:
#    unver=[nver_merge[rhi] for rhi in labels_rh_ind]#np.hstack((-100,unver1[unver1>=10242]))
#    whichhem=2
#    labels_bh_ind=copy.deepcopy(labels_rh_ind)
#labels_winhem=list()
##num_verts=np.zeros(unver[:-1].shape)
#leftright_flabels=range(2*len(unver))
#for flc in range(len(unver)):
#    print flc
#    thislabel_winner=unver[flc]#np.where(assign_vertex==unver[flc+1])[0]
#    if whichhem==1:  #left
#        thislabel_winner=thislabel_winner[thislabel_winner<10242]
#    else:
#        thislabel_winner=thislabel_winner[thislabel_winner>=10242]
#    leftright_flabels[flc*2]=thislabel_winner
#    thislabel_loser=thislabel_winner.copy()
#    for thislabii in range(len(thislabel_winner)):
#        if whichhem==1: #left
#            this_vertex=MNIlh[thislabel_winner[thislabii],:]
#            this_dist=np.sqrt(np.sum((MNIrh-this_vertex)**2,1))
#            min_dist=np.min(this_dist)
#            w_min_dist=np.argmin(this_dist)
#            while w_min_dist+10242 in thislabel_loser:#np.sum(thislabel_left==w_min_dist+10242)==1:
#                this_dist[w_min_dist]=20000;
#                min_dist=np.min(this_dist)
#                w_min_dist=np.argmin(this_dist)
#            thislabel_loser[thislabii]=w_min_dist+10242;
#        else: 
#            this_vertex=MNIrh[thislabel_winner[thislabii]-10242,:]
#            this_dist=np.sqrt(np.sum((MNIlh-this_vertex)**2,1))
#            min_dist=np.min(this_dist)
#            w_min_dist=np.argmin(this_dist)
#            while w_min_dist in thislabel_loser:#np.sum(thislabel_left==w_min_dist+10242)==1:
#                this_dist[w_min_dist]=20000;
#                min_dist=np.min(this_dist)
#                w_min_dist=np.argmin(this_dist)
#            thislabel_loser[thislabii]=w_min_dist;
#    leftright_flabels[flc*2+1]=thislabel_loser
#for flc in range(len(leftright_flabels)):
#    print flc
#    #    thislabel=unver[flc]#np.where(assign_vertex==unver[flc+1])[0]
#    #    num_verts[flc]=thislabel.shape[0]
#    thislabel_data=np.zeros((assign_vertex.shape[0],1))
#    thislabel_data[leftright_flabels[flc]]=1.0
#    thislabel_stc = mne.SourceEstimate(thislabel_data, vertices=vertices_avg, tmin=1e-3*2, tstep=1e-3*2, subject=subject_avg, verbose=None)
#    label_append=[x for x in mne.stc_to_label(thislabel_stc,src=src_avg, smooth=True,connected=False, subjects_dir=data_path) if x is not None]
#    label_append[0].smooth(smooth=2,subject='fsaverage',subjects_dir=data_path)
#    if label_append[0].hemi=='lh':#unver[flc+1]<len(assign_vertex)/2.0:
#        label_append[0].name=labels_merged[labels_bh_ind[int(flc/2)]].name[:-3]+'-lh'
#        #hemi_src = src_avg[0]
#    else:
#        label_append[0].name=labels_merged[labels_bh_ind[int(flc/2)]].name[:-3]+'-rh' 
#        #hemi_src = src_avg[1]
#    labels_winhem.append(label_append[0])
#labels=copy.deepcopy(labels_winhem)
#balanced_labels_merged=copy.deepcopy(labels)

##prepare to write split and merged labels separately
labels_merged=[label for label in labels_merged if label.name[-2:]==label.hemi]
merged_names=[label.name for label in labels_merged]
balanced_labels_merged=list()
for lcnt01 in range(len(labels_merged)):
    if merged_names[lcnt01].endswith('lh'):
        for lcnt02 in range(len(labels_merged)):
            if merged_names[lcnt02].endswith('rh') and merged_names[lcnt02][:-3]==merged_names[lcnt01][:-3]:
                balanced_labels_merged.append(labels_merged[lcnt01])
                balanced_labels_merged.append(labels_merged[lcnt02])

balanced_labels00=labels_fc+balanced_labels_merged
colors=nc(int(len(balanced_labels00)/2),bytes_=True, cmap='gist_rainbow')
colors=colors[np.random.permutation(int(len(balanced_labels00)/2))]
colors_split=colors[:int(len(labels_fc)/2),:]
fname_split='rev1_splitlabels2009_prefinal'+str(len(labels_fc))
make_parcellation_mine(src_avg,labels_fc,colors_split,subject_avg,fname_split)
## write split labels

## write merged labels
colors_merge=colors[int(len(labels_fc)/2):,:]
fname_merge='rev1_mergedlabels2009_prefinal'+str(len(balanced_labels_merged))
make_parcellation_mine(src_avg,balanced_labels_merged,colors_merge,subject_avg,fname_merge)

## write split and merged labels
fname_splitmerge='rev1_balanced_splitmerged_labels2009_prefinal'+str(len(balanced_labels00))
make_parcellation_mine(src_avg,balanced_labels00,colors,subject_avg,fname_splitmerge)


split_labels01 = mne.read_labels_from_annot(subject=subject_avg, parc=subject_avg+'_'+fname_split, subjects_dir=subjects_dir)
merged_labels01 = mne.read_labels_from_annot(subject=subject_avg, parc=subject_avg+'_'+fname_merge, subjects_dir=subjects_dir)

balanced_labels01=copy.deepcopy(balanced_labels00)#split_labels01+merged_labels01
print "stage 2 finished- add merged labels- annot file created"


#for lbl in balanced_labels01:
#    lbl.values.fill(1.0)
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
    sub_src=inverse_operator_eegmeg['src']
    fnamemf=data_path + subjects[subject_no]+'/label/'+subjects[subject_no]+'_morphed_split_aparc_2009'
    with open(fnamemf, 'rb') as f:
        labels_split=pickle.load(f)    
#    labels_split = mne.read_labels_from_annot(subject=subjects[subject_no], parc=subjects[subject_no]+'_morphed_split_aparc', subjects_dir=subjects_dir)#parc='rez_aparc_24sep15'
#    labels_split.remove(labels_split[-1])# unknown lh
#    labels_split.remove(labels_split[-1])# unknown rh
    balanced_labels_merged_copy=copy.deepcopy(balanced_labels_merged)
    morphed_merged_labels=[balanced_labels_merged_copy[jj].morph(subject_from='fsaverage', smooth=5, subject_to=subjects[subject_no], subjects_dir=data_path) for jj in range(len(balanced_labels_merged))]
    balanced_labels_subj=labels_split+morphed_merged_labels
    fname_subj='rev1_balanced_splitmerged_labels2009_prefinal'
    make_parcellation_mine(sub_src,balanced_labels_subj,colors,subjects[subject_no],fname_subj)
    labels = [label for label in balanced_labels_subj if label.name[-2:]==label.hemi]#copy.deepcopy(balanced_labels_subj)#mne.read_labels_from_annot(subject=subjects[subject_no], parc=subjects[subject_no]+'_'+fname_subj, subjects_dir=subjects_dir)#parc='rez_aparc_24sep15'
#    balanced_labels01=copy.deepcopy(labels)
#    labels = [balanced_labels01[jj].morph(subject_from='fsaverage', smooth=5, subject_to=subjects[subject_no], subjects_dir=data_path) for jj in range(len(balanced_labels01))]
    lv2=[len(label.vertices) for label in labels]
    snr = 3.0#ediv.mean()
    lambda2 = 1.0 / snr ** 2
    method = 'MNE'  # can be 'MNE' or 'sLORETA'
    mode = 'svd'
    n_svd_comp = 1
    lblname_cnt=0    
    #final ctf
    try: 
        stc_ctf_mne,expl90=psf_ctf_18Dec15.cross_talk_function(inverse_operator_eegmeg, forward, labels,method=method, lambda2=lambda2, signed=False, mode=mode, n_svd_comp=1, verbose=None)
        print "psf/ctf finished!"
        morphmat1=mne.compute_morph_matrix(subject_from=subjects[subject_no], subject_to='fsaverage', vertices_from=stc_ctf_mne.vertices, vertices_to=vertices_avg, subjects_dir=data_path)
        stc_to=stc_ctf_mne.morph_precomputed(subject_avg,vertices_avg,morphmat1, subjects[subject_no])
        if cnt==0:
            stc_all_data2=stc_to.data[:,:,np.newaxis]
        else:
            stc_all_data2=np.concatenate((stc_all_data2,stc_to.data[:,:,np.newaxis]),axis=2)
    except: 
        pass
print "stage 3 finished- ctf found for all labels- ready for the final!"

stc_mean_data2=np.mean(stc_all_data2,2)
Matx=stc_mean_data2[:,:-1]# Nv x Nlabels+1
path_mtx=data_path+'rev1_Matx_second2009.mat'
scipy.io.savemat(path_mtx,{'stc_mean_data':stc_mean_data})
Matx_zscore=np.zeros(Matx.shape)
Matx_zscore=zscore(Matx,axis=1)
assign_vertex=np.argmax(Matx_zscore, axis=1)
msort=np.sort(Matx_zscore,axis=1)
mtokeep=msort[:,-1]-msort[:,-2]
for avii in range(len(assign_vertex)):
    if Matx[avii,assign_vertex[avii]]<np.max(Matx[:,assign_vertex[avii]])/2.:
        assign_vertex[avii]=-100
assign_vertex[np.logical_or(msort[:,-1]<3, mtokeep<=1.0)]=-100
unver=np.unique(assign_vertex)
print unver
labvernum=np.zeros(unver.shape[0])
for uc in range(len(unver)):
    labvernum[uc]=np.where(assign_vertex==unver[uc])[0].shape[0]
    if labvernum[uc]<10.0:
        assign_vertex[np.where(assign_vertex==unver[uc])[0]]=-100
unver=np.unique(assign_vertex)
labvernum=np.zeros(unver.shape[0])
for uc in range(len(unver)):
    labvernum[uc]=np.where(assign_vertex==unver[uc])[0].shape[0] 
print sorted(labvernum)
print labvernum.shape

### determine winner hem
path_mtx='/imaging/rf02/parcellation/MNIlh.mat'
MNIlh_pre=scipy.io.loadmat(path_mtx,mat_dtype=True); MNIlh=np.squeeze(MNIlh_pre['MNIlh'])
path_mtx='/imaging/rf02/parcellation/MNIrh.mat'
MNIrh_pre=scipy.io.loadmat(path_mtx,mat_dtype=True); MNIrh=np.squeeze(MNIrh_pre['MNIrh'])


#labels_lh_ind=[ii1 for ii1 in range(len(unver)) if unver[ii1] in labels_lh_ind]


unver1=unver[unver>=0]
labvernum1=labvernum[1:]
labels_lh_ind=[ii1 for ii1,llcnt in enumerate(unver1) if labels[llcnt].name.endswith('lh')]
labels_rh_ind=[ii2 for ii2,rrcnt in enumerate(unver1) if labels[rrcnt].name.endswith('rh')]

if np.sum(labvernum1[labels_lh_ind])>=np.sum(labvernum1[labels_rh_ind]):#left wins
    unver=unver1[labels_lh_ind]
    whichhem=1
else:
    unver=unver1[labels_rh_ind]#np.hstack((-100,unver1[unver1>=10242]))
    whichhem=2
labels_winhem=list()
#num_verts=np.zeros(unver[:-1].shape)
leftright_flabels=range(2*len(unver))
for flc in range(len(unver)):
    print flc
    thislabel=unver[flc]#np.where(assign_vertex==unver[flc+1])[0]
    #    num_verts[flc]=thislabel.shape[0]
    thislabel_winner=np.zeros((assign_vertex.shape[0],1))
    thislabel_winner=np.where(assign_vertex==thislabel)[0]
    if whichhem==1:  #left
        thislabel_winner=thislabel_winner[thislabel_winner<10242]
    else:
        thislabel_winner=thislabel_winner[thislabel_winner>=10242]
    leftright_flabels[flc*2]=thislabel_winner
    thislabel_loser=thislabel_winner.copy()
    for thislabii in range(len(thislabel_winner)):
        if whichhem==1: #left
            this_vertex=MNIlh[thislabel_winner[thislabii],:]
            this_dist=np.sqrt(np.sum((MNIrh-this_vertex)**2,1))
            min_dist=np.min(this_dist)
            w_min_dist=np.argmin(this_dist)
            while w_min_dist+10242 in thislabel_loser:#np.sum(thislabel_left==w_min_dist+10242)==1:
                this_dist[w_min_dist]=20000;
                min_dist=np.min(this_dist)
                w_min_dist=np.argmin(this_dist)
            thislabel_loser[thislabii]=w_min_dist+10242;
        else: 
            this_vertex=MNIrh[thislabel_winner[thislabii]-10242,:]
            this_dist=np.sqrt(np.sum((MNIlh-this_vertex)**2,1))
            min_dist=np.min(this_dist)
            w_min_dist=np.argmin(this_dist)
            while w_min_dist in thislabel_loser:#np.sum(thislabel_left==w_min_dist+10242)==1:
                this_dist[w_min_dist]=20000;
                min_dist=np.min(this_dist)
                w_min_dist=np.argmin(this_dist)
            thislabel_loser[thislabii]=w_min_dist;
    leftright_flabels[flc*2+1]=thislabel_loser
print labvernum.shape
for flc in range(len(leftright_flabels)):
    print flc
    #    thislabel=unver[flc]#np.where(assign_vertex==unver[flc+1])[0]
    #    num_verts[flc]=thislabel.shape[0]
    thislabel_data=np.zeros((assign_vertex.shape[0],1))
    thislabel_data[leftright_flabels[flc]]=1.0
    thislabel_stc = mne.SourceEstimate(thislabel_data, vertices=vertices_avg, tmin=1e-3*2, tstep=1e-3*2, subject=subject_avg, verbose=None)
    label_append=[x for x in mne.stc_to_label(thislabel_stc,src=src_avg, smooth=True,connected=False, subjects_dir=data_path) if x is not None]
    label_append[0].smooth(smooth=2,subject='fsaverage',subjects_dir=data_path)
    if label_append[0].hemi=='lh':#unver[flc+1]<len(assign_vertex)/2.0:
        label_append[0].name=balanced_labels01[int(unver[flc/2])].name[:-3]+'-lh'
        #hemi_src = src_avg[0]
    else:
        label_append[0].name=balanced_labels01[int(unver[flc/2])].name[:-3]+'-rh' 
        #hemi_src = src_avg[1]
    labels_winhem.append(label_append[0])
labels=copy.deepcopy(labels_winhem)
balanced_labels=copy.deepcopy(labels)   
####################################3       
#mne.write_labels_to_annot(labels_fc, subject=subject_avg, parc=subject_avg+'_preparc_split_synccolor', subjects_dir=data_path,colormap='gist_ncar',hemi='both',overwrite=True)
#labelsf=list()
#flc1=-1
#for flc in range(len(unver)-1):
#    print flc
#    thislabel=np.where(assign_vertex==unver[flc+1])[0]
#    if thislabel.shape[0]>1.0:
#        flc1=flc1+1
#        thislabel_data=np.zeros((assign_vertex.shape[0],1))
#        thislabel_data[thislabel]=1.0
#        thislabel_stc = mne.SourceEstimate(thislabel_data, vertices=vertices_avg, tmin=1e-3*2, tstep=1e-3*2, subject=subject_avg, verbose=None)
#        label_append=[x for x in mne.stc_to_label(thislabel_stc,src=src_avg, smooth=True,connected=False, subjects_dir=data_path) if x is not None]
#        labelsf.append(label_append)
#        labelsf[flc1][0].name=balanced_labels01[unver[flc+1]].name
#        labelsf[flc1][0].color=balanced_labels01[unver[flc+1]].color
#labelsff=[label[0] for label in labelsf]
#lnames=[labelsff[lcnt01].name for lcnt01 in range(len(labelsff))]
#balanced_labels=list()
#for lcnt01 in range(len(labelsff)):
#    if lnames[lcnt01].endswith('lh'):
#        for lcnt02 in range(len(labelsff)):
#            if lnames[lcnt02].endswith('rh') and lnames[lcnt02][:-3]==lnames[lcnt01][:-3] and labelsff[lcnt02].hemi!=labelsff[lcnt01].hemi:
#                balanced_labels.append(labelsff[lcnt01])
#                balanced_labels.append(labelsff[lcnt02])
print "save the results! ^v^"
blname=[label.name for label in balanced_labels]                
#balanced_labels02 = copy.deepcopy(balanced_labels01)#mne.read_labels_from_annot(subject=subject_avg, parc='splitlabels_merged_labels_prefinal', subjects_dir=subjects_dir)
#for iii in range(len(balanced_labels02)):
#    for jjj in range(len(balanced_labels)):
#        if balanced_labels[jjj].name==balanced_labels02[iii].name:
#            balanced_labels[jjj].color=balanced_labels02[iii].color

#mne.write_labels_to_annot(balanced_labels, subject=subject_avg, parc=subject_avg+'_finalpostparc_smooth_split_merge', subjects_dir=data_path,colormap='gist_ncar',hemi='both',overwrite=True)
balanced_labels00c=balanced_labels00[::2]
clrs=np.asarray([colors[jj,:] for jj in range(len(balanced_labels00c)) if balanced_labels00c[jj].name in blname])
#balanced_labels00 = mne.read_labels_from_annot(subject=subject_avg, parc='balanced_splitmerged_labels_prefinal', subjects_dir=subjects_dir)
balanced_labels00_final=[label for label in balanced_labels00 if label.name in blname]
split_labels01name=[label.name for label in split_labels01]
merge_labels01name=[label.name for label in merged_labels01]

balanced_labelssplit_final=[label for label in balanced_labels00_final if label.name in split_labels01name]
balanced_labelsmerge_final=[label for label in balanced_labels00_final if label.name in merge_labels01name]

fname_blm='rev1_finalpreparc2009_smooth_stddiff_mergeonly_halfmax'+str(len(balanced_labelsmerge_final))
make_parcellation_mine(src_avg,balanced_labelsmerge_final,clrs[len(balanced_labelssplit_final)/2:,:],subject_avg,fname_blm)

fname_bls='rev1_finalpreparc2009_smooth_stddiff_splitonly_halfmax'+str(len(balanced_labelssplit_final))
make_parcellation_mine(src_avg,balanced_labelssplit_final,clrs[:len(balanced_labelssplit_final)/2,:],subject_avg,fname_bls)

fname_final='rev1_finalpostparc2009_smooth_stddiff_split_merge_halfmax_'+str(len(balanced_labels))
make_parcellation_mine(src_avg,balanced_labels,clrs,subject_avg,fname_final)
lv2=[len(label.vertices) for label in balanced_labels]    
col_fname='rev1_colour_neuroanat_winnerhem_'+str(len(balanced_labels))+'_2009.mat'
scipy.io.savemat(col_fname,{'clrs':clrs})
fname_blb='rev1_finalpreparc2009_smooth_stddiff_splitmerge_halfmax'
make_parcellation_mine(src_avg,balanced_labels00_final,clrs,subject_avg,fname_blb)

#mne.write_labels_to_annot(balanced_labels00_final, subject=subject_avg, parc=subject_avg+'_finalpreparc_smooth_split_merge', subjects_dir=data_path,colormap='gist_ncar',hemi='both',overwrite=True)


#path_c=data_path+'colormatrix.mat'
#scipy.io.savemat(path_c,{'rs':rs})
#lfname=[labelsff[ii].name for ii in range(len(labelsff))]
#lfname_path=data_path+subject_from+'/labels_modified_names_split_merge.mat'
#scipy.io.savemat(lfname_path,{'lfname':lfname})

