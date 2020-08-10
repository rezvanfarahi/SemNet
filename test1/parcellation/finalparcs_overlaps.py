# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 14:01:27 2016

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
reload(mne)
from mne.minimum_norm import (read_inverse_operator, psf_ctf_28Oct15)
import matplotlib.pyplot as plt
import scipy.io
from mne.label import _n_colors as nc
from mne.label import _read_annot as mnera
from mne.label import _write_annot as mnewa
import copy
#data_path the main general path, e.g. /imaging/rf02/TypLexMEG 
#list_all subdirectory under the general path where all the preprocessed data files including fwd and inverse reside e.g. data_path/meg10_0378/101209/
#subjects: subdirectory under the general path where all the preprocessed MRI files, BEM surface, etc. reside e.g. data_path/SS16_Meg11_0104/
def parcellation_leakage_mat(data_path, list_all, subjects, fwd_name, inv_name, labelss, subject_avg, vertices_avg):
    ii=-1
    for cnt, meg in enumerate(list_all): 
        print cnt
        fname_fwd = data_path + meg + fwd_name#'ico5_forward_5-3L-EMEG-fwd.fif'
        fname_inv = data_path + meg + inv_name#'InvOp_ico5newreg_fft_1_48_clean_ica_EMEG-inv.fif'
        vertices_to = [np.arange(10242), np.arange(10242)]
        subject_no=cnt
    
        forward = mne.read_forward_solution(fname_fwd,surf_ori=True, force_fixed=True)#, surf_ori=True)
    
        # read inverse operators
        inverse_operator_eegmeg = read_inverse_operator(fname_inv)
        #morph fsaverage labels to each individual space to compute ctf
        labels = [labelss[jj].morph(subject_from='fsaverage', smooth=5, subject_to=subjects[subject_no], subjects_dir=data_path) for jj in range(len(labelss))]
    
        snr = 3.0#ediv.mean()
        lambda2 = 1.0 / snr ** 2
        method = 'MNE'  # can be 'MNE' or 'sLORETA'
        mode = 'svd'
        n_svd_comp = 1                
        subject_from = subjects[subject_no]
        subject_to = 'fsaverage'
        # find ctf using svd approach, if didn't converge for a subject skip the subject
        try:
            stc_ctf_mne=psf_ctf_28Oct15.cross_talk_function(inverse_operator_eegmeg, forward, labels=labels,method=method, lambda2=lambda2, signed=False, mode=mode, n_svd_comp=n_svd_comp, verbose=None)
            ii=ii+1
            morphmat1=mne.compute_morph_matrix(subject_from, subject_to, stc_ctf_mne.vertices, vertices_to, subjects_dir=data_path)
            stc_to=stc_ctf_mne.morph_precomputed(subject_to,vertices_to,morphmat1, subject_from)
            if ii==0:
                stc_all_data=stc_to.data[:,:,np.newaxis]
            else:
                stc_all_data=np.concatenate((stc_all_data,stc_to.data[:,:,np.newaxis]),axis=2)
        except:
            pass#stc_ctf_mne=psf_ctf_28Oct15.cross_talk_function(inverse_operator_eegmeg, forward, labels=labels,method=method, lambda2=lambda2, signed=False, mode='mean', verbose=None)
        
        #morph ctf back to the average space to take the average over subjects
        
    stc_mean_data=np.mean(stc_all_data,axis=2)
    all_stc = mne.SourceEstimate(stc_mean_data[:,:-1]/stc_mean_data[:,-1][:,np.newaxis], vertices=vertices_avg, tmin=1e-3*2, tstep=1e-3*2, subject=subject_avg, verbose=None)
    ctf_to=np.zeros((len(labelss),len(labelss)))
    
    for ii in range(len(labelss)):
        ctf_to[ii,:]=np.mean(all_stc.in_label(labelss[ii]).data,axis=0) #each row: cross talk received at each label
    return ctf_to
    
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
##################!!!!NOTE names RG and anat are used reversed
data_path = '/imaging/rf02/TypLexMEG/'
parc_path=data_path+'parcellation/final_results/neuroanat_merge_nomirror/'#neuroanat_merge_nomirror/'
out_path = '/imaging/rf02/'
RG_name='fsaverage_finalpostparc_smooth_stddiff_split_merge_halfmax_86_badremoved_74'#'fsaverage_finalpostparc_smooth_stddiff_split_merge_halfmax_86'#'fsaverage_RG_finalpostparc_smooth_stddiff_split_merge_halfmax_82'#'fsaverage_finalpostparc_smooth_stddiff_split_merge_halfmax_86'
anat_name='fsaverage_finalpostparc2009_smooth_stddiff_split_merge_halfmax_98_badremoved_74'#'fsaverage_sorted_RG_finalpostparc_smooth_stddiff_split_merge_halfmax_82_badremoved_70'
labels_anat1 = mne.read_labels_from_annot(subject='fsaverage', parc=anat_name, subjects_dir=data_path)#parc='rez_aparc_24sep15'
labels_anat=labels_anat1[::2]

for lbl in labels_anat:
    lbl.values.fill(1.0)

labels_RG1 = mne.read_labels_from_annot(subject='fsaverage', parc=RG_name, subjects_dir=data_path)#parc='rez_aparc_24sep15'
labels_RG=labels_RG1[::2]

for lbl in labels_RG:
    lbl.values.fill(1.0)
stc_file=data_path + 'postcentral_CTF'
stc_obj=mne.read_source_estimate(stc_file)
anat_nvert=np.asarray([len(label.vertices) for label in labels_anat])
anat_vert_sort=np.argsort(anat_nvert)[-1::-1]
not_avail=list()
assign_list=range(len(labels_anat))
overlap_mat=np.zeros((len(labels_anat),len(labels_RG)))
for anatcnt,anatvert in enumerate(anat_vert_sort):
    thisanat_verts=labels_anat[anatvert].vertices#stc_obj.in_label(labels_anat[anatvert]).vertices[0]
    this_overlap=np.ones((len(labels_RG),1))*-1
    for RGcnt, thisRG in enumerate(labels_RG):        
        if RGcnt not in not_avail:
            thisRG_verts=thisRG.vertices#stc_obj.in_label(thisRG).vertices[0]
            this_overlap[RGcnt]=len(np.intersect1d(thisanat_verts,thisRG_verts))
        overlap_mat[anatvert,RGcnt]=len(np.intersect1d(thisanat_verts,thisRG_verts))/np.float(len(thisanat_verts))
    if anatcnt<len(labels_RG):
        assign_list[anatvert]=np.argmax(this_overlap)
        not_avail.append(np.argmax(this_overlap))
    else:
        assign_list[anatvert]=-1
not_avail2=np.array([])
not_avail3=-1*np.ones((len(anat_vert_sort),))

for anatcnt,anatvert in enumerate(anat_vert_sort):#ii in range(overlap_mat.shape[0]):
    thisrow=overlap_mat[anatvert,:]
    thismax=np.argsort(thisrow)[-1::-1]
    thisov=np.in1d(thismax,not_avail2)
    not_avail2=np.append(not_avail2,thismax[~thisov][0])
    not_avail3[anatvert]=thismax[~thisov][0]
long_len=np.arange(overlap_mat.shape[1])
still_avail=~np.in1d(long_len,not_avail3.astype(int))
not_avail3=np.append(not_avail3,long_len[still_avail])
order_match=np.zeros((2*not_avail3.shape[0],))
order_match[::2]=not_avail3*2
order_match[1::2]=not_avail3*2+1
order_match=order_match.astype(int)
path_to=out_path+'parcellation/leakage_pattern/anat_anat2009_sameorder.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot
scipy.io.savemat(path_to,{'order_match':order_match})

plt.imshow(overlap_mat[:,not_avail3.astype(int)], cmap='hot', aspect='equal', interpolation="nearest")#, vmin=0, vmax=1)
plt.xlabel('Split and Merge labels')
plt.ylabel('Split and Merge 2009 labels')
plt.title('Percentage of Label Overlaps')

cbar=plt.colorbar()
#plt.show()
#cbar.ax.tick_params(labelsize=40)
#labelss_name=[str(label.name)+"\n" for label in labels_final]
#s=" "
#plt.ylabel(s.join(labelss_name),rotation=0, horizontalalignment='center', verticalalignment='center',fontsize=4, labelpad=10)
#plt.xlabel(s.join(labelss_name),rotation=90,fontsize=6, labelpad=40)
#plt.subplots_adjust(left=0.1, right=0.9, top=0.95, bottom=0.4)
#out_name=out_path+'parcellation/leakage_pattern/'+new_fname+'.jpg'
##out_name2=out_path+'parcellation/leakage_pattern/'+'labels.tif'
#plt.savefig(out_name, figsize=(100,100), facecolor='white',format='jpg',dpi=300)#,tight_layout=True)
#plt.close()
#
new_colors=np.zeros((len(not_avail),4))
new_colors[:,3]=1.0
for lfcnt,thislabelc in enumerate(not_avail):#lfcnt,thislabelf in enumerate(labels_RG):
    this_color=np.asarray(labels_RG[thislabelc].color)[:3]
    new_colors[anat_vert_sort[lfcnt],:3]=this_color
final_colors=(new_colors*255).astype('uint8')
final_name='finalpostparc2009_smooth_stddiff_split_merge_halfmax_98_badremoved_74_cmatch'
fwd_name='ico5_forward_5-3L-EMEG-fwd.fif'
inv_name='InvOp_ico5newreg_fft_1_48_clean_ica_EMEG-inv.fif'
vertices_avg = [np.arange(10242), np.arange(10242)]
subject_avg='fsaverage'
srcin='/imaging/rf02/TypLexMEG/fsaverage/bem/fsaverage_dist-ico-5-src.fif'
src_avg = mne.read_source_spaces(srcin)
make_parcellation_mine(src_avg,labels_anat1,final_colors,subject_avg,final_name)
