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

labelss = mne.read_labels_from_annot(subject='fsaverage', parc='aparc', subjects_dir=subjects_dir)#parc='rez_aparc_24sep15'

#labelss.remove(labelss[-3])
#labelss.remove(labelss[-2])
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

label_colors = [label.color for label in labelss]
# read label(s)
for lbl in labelss:
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
    labels = mne.read_labels_from_annot(subject=subjects[subject_no], parc='morphed_aparc_21Dec15', subjects_dir=subjects_dir)#parc='rez_aparc_24sep15'
    labelname=[labels[ii].name for ii in range(len(labels))]
    labels.remove(labels[labelname.index('frontalpole-lh')])
    labels.remove(labels[labelname.index('frontalpole-lh')])#don't panic! it's actually right hemishphere
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
    path_c=data_path+'explained80all.mat'
    ltc=scipy.io.loadmat(path_c,mat_dtype=True); explained80=np.squeeze(ltc['EXP80'])
    explained80=np.delete(explained80,[labelname.index('frontalpole-lh'),labelname.index('frontalpole-lh')+1],0)
    llnew=[labels[ii] for ii in range(len(labels)) if explained80[ii]==1 or len(labels[ii].vertices)<10]
    lln1=[labels[ii].split(explained80[ii], subjects_dir=data_path)[:] for ii in range(len(labels)) if explained80[ii]>1 and len(labels[ii].vertices)>=10]
    llnn=list()
    for ii in range(len(lln1)):
        llnn+=lln1[ii]
    
    labels_backup=copy.deepcopy(labels)
    labels=llnn+llnew
    label_names=[labels[lcn3].name for lcn3 in range(len(labels))]
    ordered_indices=np.zeros((len(label_names),1))
    lcntl=-1
    for lcnt3 in range(len(label_names)):
        if label_names[lcnt3].endswith('lh'):
            lcntl=lcntl+1
            thislabelnameR=label_names.index(label_names[lcnt3][:-3]+'-rh')
            this_indices=np.array([lcnt3,thislabelnameR])
            ordered_indices[2*lcntl]=lcnt3; ordered_indices[2*lcntl+1]=thislabelnameR
    labels=[labels[jjk] for jjk in np.squeeze(ordered_indices.astype(int))]
    subject_from = subjects[subject_no]
    subject_to = 'fsaverage'
    colors=nc(len(labels),bytes_=False, cmap='gist_ncar')#Note! left right different colors, not important here.
    colors=colors[np.random.permutation(int(len(labels)))]
    for lcnt1 in range(len(labels)):
        labels[lcnt1].color=colors[lcnt1,:]
    mne.write_labels_to_annot(labels, subject=subject_from, parc=subject_from+'_morphed_split_aparc', subjects_dir=data_path,colormap='gist_ncar',hemi='both',overwrite=True)
    labels = mne.read_labels_from_annot(subject=subjects[subject_no], parc=subjects[subject_no]+'_morphed_split_aparc', subjects_dir=subjects_dir)#parc='rez_aparc_24sep15'
    labels.remove(labels[-1])# unknown rh
    labels.remove(labels[-1])# unknown lh    
    stc_ctf_mne, ep80=psf_ctf_18Dec15.cross_talk_function(inverse_operator_eegmeg, forward, labels,method=method, lambda2=lambda2, signed=False, mode=mode, n_svd_comp=1, verbose=None)
    lv=[len(labels[vc].vertices) for vc in range(len(labels))]
    print "psf/ctf finished!"
    #mne.write_labels_to_annot(labels, subject=subject_from, parc=subject_from+'_preparc_smooth_split', subjects_dir=data_path,colormap='gist_ncar',hemi='both',overwrite=True)
    #morphmat1=mne.compute_morph_matrix(subject_from, subject_to, stc_ctf_mne.vertices, vertices_to, subjects_dir=data_path)
    stc_to=stc_ctf_mne.copy()#stc_ctf_mne.morph_precomputed(subject_to,vertices_to,morphmat1, subject_from)
    Matx=stc_to.data[:,:-1]# Nv x Nlabels+1
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
    num_oldv=np.zeros((len(uv)-1,1))
    for olcnt in range(len(uv)-1):
        num_oldv[olcnt]=np.where(assign_vertex==uv[olcnt+1])[0].shape[0]
    
    merge_candidates=np.where(np.logical_and(np.logical_and(msort[:,-1]>=3,msort[:,-2]>=3), mtokeep<=1))
    margsort=np.argsort(Matx_zscore,axis=1)
    merge_labels1=margsort[merge_candidates[0],-2:]
    merge_labels=np.sort(merge_labels1,axis=1)
    new_labels=[tuple(row) for row in merge_labels]
    newlabels=np.unique(new_labels)
    lablen=len(labels)
    nlsF=newlabels[newlabels[:,0]%2==newlabels[:,1]%2,:]
    merge_keep=np.zeros(nlsF.shape)
    mcnt=-1
    for lcnt in range(nlsF.shape[0]):
        for rcnt in range(nlsF.shape[0]):
            if nlsF[lcnt,0]%2==0: ##lh labels
                if np.array_equal(nlsF[lcnt,:],nlsF[rcnt,:]-1):
                    mcnt=mcnt+1
                    merge_keep[mcnt,:]=nlsF[lcnt,:]
                
    merge_keepL=merge_keep[:mcnt+1,:]; merge_keepR=merge_keepL+ 1        
    merge_keepB_new=np.concatenate((merge_keepL,merge_keepR),axis=0)
    if cnt==0:
        merge_keepB_old=merge_keepB_new
            
    merge_keep_all=np.zeros(merge_keepB_old.shape)
    mcnt=-1
    for lcnt in range(merge_keepB_old.shape[0]):
        for rcnt in range(merge_keepB_new.shape[0]):
                if np.array_equal(merge_keepB_old[lcnt,:],merge_keepB_new[rcnt,:]):
                    mcnt=mcnt+1
                    merge_keep_all[mcnt,:]=merge_keepB_old[lcnt,:]
    merge_keepB_old=merge_keep_all[:mcnt+1,:]
    print merge_keepB_old.shape
merge_final=np.zeros(merge_keepB_old.shape)
mcnt=-1
for lcnt in range(merge_keepB_old.shape[0]):
    for rcnt in range(merge_keepB_old.shape[0]):
        if merge_keepB_old[lcnt,0]%2==0: ##lh labels
            if np.array_equal(merge_keepB_old[lcnt,:],merge_keepB_old[rcnt,:]-1):
                mcnt=mcnt+1
                merge_final[mcnt,:]=merge_keepB_old[lcnt,:]
keepL=merge_final[:mcnt+1,:];   keepR=keepL+ 1 
newlabels=np.concatenate((keepL,keepR),axis=0) 
path_c=data_path+'labelpairs_merge80_allsubjs2.mat'
scipy.io.savemat(path_c,{'newlabels':newlabels})   
#    labels_merged=list()
#    lblname_cnt=-1
#    for labcnt in range(newlabels.shape[0]):
#        print labcnt
#        thislabel=np.where(np.logical_and(merge_labels[:,0]==newlabels[labcnt,0],merge_labels[:,1]==newlabels[labcnt,1]))[0]
#        if thislabel.shape[0]>=25:
#            lblname_cnt=lblname_cnt+1
#            thislabel_data=np.zeros((assign_vertex.shape[0],1))
#            thislabel_data[thislabel]=1.0
#            thislabel_stc = mne.SourceEstimate(thislabel_data, vertices=stc_ctf_mne.vertices, tmin=1e-3*2, tstep=1e-3*2, subject=subject_from, verbose=None)
#            label_append=[x for x in mne.stc_to_label(thislabel_stc,src=sub_src, smooth=True,connected=False, subjects_dir=data_path) if x is not None]
#            labels_merged.append(label_append)
#            labels_merged[lblname_cnt][0].name=labels[newlabels[labcnt,1]].name[:-3]+'_'+labels[newlabels[labcnt,0]].name
#    unver=np.unique(assign_vertex)
#    
#    labvernum=np.zeros(unver.shape[0])
#    for uc in range(len(unver)):
#        labvernum[uc]=np.where(assign_vertex==unver[uc])[0].shape[0]
#        if labvernum[uc]<1.0:
#            assign_vertex[np.where(assign_vertex==unver[uc])[0]]=-100
#    unver=np.unique(assign_vertex)
#    print unver
#    labvernum=np.zeros(unver.shape[0])
#    for uc in range(len(unver)):
#        labvernum[uc]=np.where(assign_vertex==unver[uc])[0].shape[0]        
#    print sorted(labvernum)
#    print labvernum.shape
#    for lcnt in range(len(unver)-1):
#        thislabel=np.where(assign_vertex==unver[lcnt+1])[0]
#        thislabel_data=np.zeros((assign_vertex.shape[0],1))
#        thislabel_data[thislabel]=1.0
#        thislabel_stc = mne.SourceEstimate(thislabel_data, vertices=stc_ctf_mne.vertices, tmin=1e-3*2, tstep=1e-3*2, subject=subject_from, verbose=None)
#        label_append=[x for x in mne.stc_to_label(thislabel_stc,src=sub_src, smooth=True,connected=False, subjects_dir=data_path) if x is not None]
#        labels_merged.append(label_append)
#        labels_merged[lblname_cnt+lcnt+1][0].name=labels[unver[lcnt+1]].name
#    labels=[label[0] for label in labels_merged]
#    
#    
#    stc_ctf_mne,expl80=psf_ctf_18Dec15.cross_talk_function(inverse_operator_eegmeg, forward, labels,method=method, lambda2=lambda2, signed=False, mode=mode, n_svd_comp=1, verbose=None)
#    print "psf/ctf finished!"
#    stc_to=stc_ctf_mne.copy()#stc_ctf_mne.morph_precomputed(subject_to,vertices_to,morphmat1, subject_from)
#    Matx=stc_to.data[:,:-1]# Nv x Nlabels+1
#    Matx_zscore=np.zeros(Matx.shape)
#    Matx_zscore=zscore(Matx,axis=1)
#    assign_vertex=np.argmax(Matx_zscore, axis=1)
#    msort=np.sort(Matx_zscore,axis=1)
#    mtokeep=msort[:,-1]-msort[:,-2]
#    for avii in range(len(assign_vertex)):
#        if Matx[avii,assign_vertex[avii]]<np.max(Matx[:,assign_vertex[avii]])/4.:
#            assign_vertex[avii]=-100
#    assign_vertex[np.logical_or(msort[:,-1]<3, mtokeep<=1)]=-100
#    unver=np.unique(assign_vertex)
#    print unver
#    labvernum=np.zeros(unver.shape[0])
#    for uc in range(len(unver)):
#        labvernum[uc]=np.where(assign_vertex==unver[uc])[0].shape[0] 
#    print sorted(labvernum)
#    print labvernum.shape
#    labels_final=[labels[ii] for ii in unver[1:]]
#    label_names=[label.name for label in labels_final]
#    
#    rs=nc(len(labels_final),bytes_=False, cmap='gist_ncar')
#    rs=rs[np.random.permutation(len(labels_final))]
#    import copy
#    labels_removeov=copy.deepcopy((labels))
#    for cnt1 in range(len(unver[1:])):
#        print cnt1
#        for cnt2 in range(len(unver[1:])):
#            if cnt1<cnt2:
#                label1ver=labels_removeov[unver[cnt1+1]].vertices
#                label2ver=labels_removeov[unver[cnt2+1]].vertices
#                vnoverlap=np.in1d(label2ver,label1ver,invert=True)#voverlap=np.intersect1d(label1ver,label2ver)
#                #overlap_index=np.where(voverlap)[0]
#                label2ver_new=label2ver[vnoverlap]#np.delete(label2ver,overlap_index)
#                jjh=labels_removeov[unver[cnt2+1]]
#                jjh.vertices=label2ver_new
#        labels_removeov[unver[cnt1+1]].color=rs[cnt1]
#    labels_semifinal=[labels_removeov[ii] for ii in unver[1:]]
#    mne.write_labels_to_annot(labels_semifinal, subject=subject_from, parc=subject_from+'_test_postparc_split_merge', subjects_dir=data_path,colormap='gist_ncar',hemi='both',overwrite=True)
#    for cc1 in range(len(labelss)):
#        for cc2 in range(len(labels_semifinal)):
#            if labelss[cc1].name==labels_semifinal[cc2].name:
#                labelss[cc1].color=labels_semifinal[cc2].color
#    mne.write_labels_to_annot(labelss, subject='fsaverage', parc='test_fsaverage_preparc_split_merge', subjects_dir=data_path,colormap='gist_ncar',hemi='both',overwrite=True)
#    labelsf=list()
#    for flc in range(len(unver)-1):
#        print flc
#        thislabel=np.where(assign_vertex==unver[flc+1])[0]
#        thislabel_data=np.zeros((assign_vertex.shape[0],1))
#        thislabel_data[thislabel]=1.0
#        thislabel_stc = mne.SourceEstimate(thislabel_data, vertices=stc_ctf_mne.vertices, tmin=1e-3*2, tstep=1e-3*2, subject=subject_from, verbose=None)
#        label_append=[x for x in mne.stc_to_label(thislabel_stc,src=sub_src, smooth=True,connected=False, subjects_dir=data_path) if x is not None]
#        labelsf.append(label_append)
#        labelsf[flc][0].name=labels_semifinal[flc].name
#        labelsf[flc][0].color=labels_semifinal[flc].color
#    labelsff=[label[0] for label in labelsf]
#    mne.write_labels_to_annot(labelsff, subject=subject_from, parc=subject_from+'_test_finalpostparc_smooth_split_merge', subjects_dir=data_path,colormap='gist_ncar',hemi='both',overwrite=True)
##path_c=data_path+'colormatrix.mat'
##scipy.io.savemat(path_c,{'rs':rs})
##lfname=[labelsff[ii].name for ii in range(len(labelsff))]
##lfname_path=data_path+subject_from+'/labels_modified_names_split_merge.mat'
##scipy.io.savemat(lfname_path,{'lfname':lfname})
