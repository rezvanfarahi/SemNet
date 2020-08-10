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
from mne.minimum_norm import (read_inverse_operator, make_inverse_operator, write_inverse_operator, apply_inverse, apply_inverse_epochs,psf_ctf_new, point_spread_function,psf_ctf_28Oct15)
from mne.minimum_norm.inverse import (prepare_inverse_operator)
import scipy.io
from mne import find_events
from mne.connectivity import seed_target_indices, spectral_connectivity
from mne import read_labels_from_annot
from mne.label import _n_colors as nc
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

subject_inds=[0]#,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
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
labelss = mne.read_labels_from_annot(subject='fsaverage', parc='rez_aparc_24sep15', subjects_dir=subjects_dir)
#labelss.remove(labelss[-3])
#labelss.remove(labelss[-2])
labelss.remove(labelss[-1])
print len(labelss)
llk=labelss[0]+labelss[-2]
lrk=labelss[1]+labelss[-1]
llk.name=labelss[0].name
lrk.name=labelss[1].name
labelss.remove(labelss[-1])
labelss.remove(labelss[-1])
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
    print cnt, "preload variables"
    subjects_dir = data_path 
    fname_fwd = data_path + meg + 'ico5_forward_5-3L-EMEG-fwd.fif'
    fname_inv = inv_path + meg + 'InvOp_ico5newreg_fft_1_48_clean_ica_EMEG-inv.fif'
    vertices_to = [np.arange(10242), np.arange(10242)]
    subject_no=subject_inds[cnt]

    forward = mne.read_forward_solution(fname_fwd,surf_ori=True, force_fixed=True)#, surf_ori=True)

    # read inverse operators
    inverse_operator_eegmeg = read_inverse_operator(fname_inv)
    
    snr = 3.0#ediv.mean()
    lambda2 = 1.0 / snr ** 2
    method = 'MNE'  # can be 'MNE' or 'sLORETA'
    mode = 'svd'
    n_svd_comp = 1
    lblname_cnt=0
    sub_src=inverse_operator_eegmeg['src']
    # compute ctf for all vertices
    stc_ctf_mne=psf_ctf_28Oct15.cross_talk_function(inverse_operator_eegmeg, forward, labels=None,method=method, lambda2=lambda2, signed=False, mode=mode, n_svd_comp=1, verbose=None)
    print "psf/ctf finished!"

    subject_from = subjects[subject_no]
    subject_to = 'fsaverage'
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

    unver=unique(assign_vertex)
#    vcolor=nc(len(unver[:-1]),bytes_=False, cmap='gist_ncar')
    labelsv=list()
    num_verts=np.zeros(unver[:-1].shape)
    for flc in range(len(unver)-1):
        print flc
        thislabel=np.where(assign_vertex==unver[flc+1])[0]
        num_verts[flc]=thislabel.shape[0]
        thislabel_data=np.zeros((assign_vertex.shape[0],1))
        thislabel_data[thislabel]=1.0
        thislabel_stc = mne.SourceEstimate(thislabel_data, vertices=stc_ctf_mne.vertices, tmin=1e-3*2, tstep=1e-3*2, subject=subject_from, verbose=None)
        label_append=[x for x in mne.stc_to_label(thislabel_stc,src=sub_src, smooth=True,connected=False, subjects_dir=data_path) if x is not None]
        labelsv.append(label_append)
        if unver[flc+1]<len(assign_vertex)/2.0:
            labelsv[flc][0].name='L'+str(unver[flc+1])+'-lh'
        else:
            labelsv[flc][0].name='L'+str(int(unver[flc+1]-len(assign_vertex)/2))+'-rh'
        #labelsv[flc][0].color=vcolor[flc]    
    labelsvv=[label[0] for label in labelsv]
    mean_vert=num_verts.mean()
    mean_vert_pre=1.0
    mean_vert_change=(mean_vert-mean_vert_pre)/mean_vert
    labels=copy.deepcopy(labelsvv)
    whilec=-1
    while mean_vert_change>0.05:
        whilec=whilec+1
        print whilec
        stc_ctf_mne=psf_ctf_28Oct15.cross_talk_function(inverse_operator_eegmeg, forward, labels=labels,method=method, lambda2=lambda2, signed=False, mode=mode, n_svd_comp=1, verbose=None)
        print "psf/ctf finished!"
        stc_to=stc_ctf_mne.copy()#stc_ctf_mne.morph_precomputed(subject_to,vertices_to,morphmat1, subject_from)
        Matx=stc_to.data[:,:-1]# Nv x Nlabels+1
        Matx_zscore=np.zeros(Matx.shape)
        Matx_zscore=zscore(Matx,axis=1)
        assign_vertex=np.argmax(Matx_zscore, axis=1)
        msort=np.sort(Matx_zscore,axis=1)
        mtokeep=msort[:,-1]-msort[:,-2]
        for avii in range(len(assign_vertex)):
            if Matx[avii,assign_vertex[avii]]<np.max(Matx[:,assign_vertex[avii]])/4.:
                assign_vertex[avii]=-100
        
        assign_vertex[logical_or(msort[:,-1]<3, mtokeep<=1)]=-100
        merge_candidates=np.where(logical_and(logical_and(msort[:,-1]>=3,msort[:,-2]>=3), mtokeep<=1))
        margsort=np.argsort(Matx_zscore,axis=1)
        merge_labels1=margsort[merge_candidates[0],-2:]
        merge_labels=np.sort(merge_labels1,axis=1)
        new_labels=[tuple(row) for row in merge_labels]
        newlabels=np.unique(new_labels)
        labels_merged=list()
        lblname_cnt=-1
        ## for new merged labels
        for labcnt in range(newlabels.shape[0]):
            print labcnt
            thislabel=np.where(logical_and(merge_labels[:,0]==newlabels[labcnt,0],merge_labels[:,1]==newlabels[labcnt,1]))[0]
            if thislabel.shape[0]>=mean_vert:
                lblname_cnt=lblname_cnt+1
                thislabel_data=np.zeros((assign_vertex.shape[0],1))
                thislabel_data[thislabel]=1.0
                thislabel_stc = mne.SourceEstimate(thislabel_data, vertices=stc_ctf_mne.vertices, tmin=1e-3*2, tstep=1e-3*2, subject=subject_from, verbose=None)
                label_append=[x for x in mne.stc_to_label(thislabel_stc,src=sub_src, smooth=True,connected=False, subjects_dir=data_path) if x is not None]
                labels_merged.append(label_append)
                labels_merged[lblname_cnt][0].name=labels[newlabels[labcnt,1]].name[:-3]+'_'+labels[newlabels[labcnt,0]].name
        ## for the labels we already had
        unver=unique(assign_vertex)        
        labvernum=np.zeros(unver.shape[0])
        for uc in range(len(unver)):
            labvernum[uc]=np.where(assign_vertex==unver[uc])[0].shape[0]
            if labvernum[uc]<1.0:
                assign_vertex[np.where(assign_vertex==unver[uc])[0]]=-100
        unver=unique(assign_vertex)
        #print unver
        labvernum=np.zeros(unver.shape[0])
        for uc in range(len(unver)):
            labvernum[uc]=np.where(assign_vertex==unver[uc])[0].shape[0]        
        #print sorted(labvernum)
        print labvernum.shape
        for lcnt in range(len(unver)-1):
            thislabel=np.where(assign_vertex==unver[lcnt+1])[0]
            thislabel_data=np.zeros((assign_vertex.shape[0],1))
            thislabel_data[thislabel]=1.0
            thislabel_stc = mne.SourceEstimate(thislabel_data, vertices=stc_ctf_mne.vertices, tmin=1e-3*2, tstep=1e-3*2, subject=subject_from, verbose=None)
            label_append=[x for x in mne.stc_to_label(thislabel_stc,src=sub_src, smooth=True,connected=False, subjects_dir=data_path) if x is not None]
            labels_merged.append(label_append)
            labels_merged[lblname_cnt+lcnt+1][0].name=labels[unver[lcnt+1]].name
        labels=[label[0] for label in labels_merged]
        num_verts=np.zeros((len(labels),1))
        for labc in range(len(labels)):
            num_verts[labc]=labels[labc].vertices.shape[0]
        mean_vert_pre=mean_vert.copy()
        mean_vert=num_verts.mean()
        mean_vert_change=(mean_vert-mean_vert_pre)/mean_vert
        print len(labels)
    ### end of the while loop
    rs=nc(len(labels),bytes_=False, cmap='gist_ncar')
    rs=rs[np.random.permutation(len(labels))]
    for colorc in range(len(labels)):
        labels[colorc].color=rs[colorc]
    mne.write_labels_to_annot(labels, subject=subject_from, parc=subject_from+'_vertices2', subjects_dir=data_path,colormap='gist_ncar',hemi='both',overwrite=True)
    
#    labelsff=[label[0] for label in labelsf]
#    mne.write_labels_to_annot(labelsff, subject=subject_from, parc=subject_from+'_finalpostparc_nonsmooth', subjects_dir=data_path,colormap='gist_ncar',hemi='both',overwrite=True)
#    path_c=data_path+'colormatrix.mat'
#    scipy.io.savemat(path_c,{'rs':rs})
#    lfname=[labelsff[ii].name for ii in range(len(labelsff))]
#    lfname_path=data_path+subject_from+'/labels_modified_names.mat'
#    scipy.io.savemat(lfname_path,{'lfname':lfname})
    
