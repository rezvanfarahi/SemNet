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
reload(mne)
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

def make_parcellation_mine_vert(src,labels,colors,subject,fname):
    annotl=np.zeros((src[0]['np'],1))
    annotr=np.zeros((src[1]['np'],1))
    lh_labels=[label for label in labels if label.hemi=='lh']
    rh_labels=[label for label in labels if label.hemi=='rh']
    CTAB_L=np.zeros((len(lh_labels),5))
    CTAB_R=np.zeros((len(rh_labels),5))
    colors_lh=colors[:len(lh_labels),:]
    colors_rh=colors[len(lh_labels):,:]
    annot_id_coding = np.array((1, 2 ** 8, 2 ** 16))
    annot_id_lh=np.sum(colors_lh[:,:3] * annot_id_coding, axis=1)
    annot_id_rh=np.sum(colors_rh[:,:3] * annot_id_coding, axis=1)
    lhcnt=-1; rhcnt=-1
    for lcntn in range(len(labels)):
        print lcntn
        if labels[lcntn].name.endswith('lh'):
            lhcnt=lhcnt+1
            annotl[labels[lcntn].vertices]=annot_id_lh[lhcnt]
        else:
            rhcnt=rhcnt+1
            annotr[labels[lcntn].vertices]=annot_id_rh[rhcnt]
    CTAB_L[:,:4]=colors_lh; CTAB_L[:,4]=annot_id_lh
    CTAB_R[:,:4]=colors_rh; CTAB_R[:,4]=annot_id_rh
    #label_names=[labels[lh].name[:-3] for lh in range(len(labels))]
    NAMES_L=[label.name[:-3] for label in labels if label.hemi=='lh']#label_names[::2]
    NAMES_R=[label.name[:-3] for label in labels if label.hemi=='rh']#label_names[1::2]
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

ll = []
for ss in subject_inds:
 ll.append(list_all[ss])


print "ll:"
print ll
ii=-1

stc_all=range(17)
assign_vertex=np.zeros((20484,17))


#for cnt, meg in enumerate(ll):
#    ii=ii+1
#    print cnt, "preload variables"
#    subjects_dir = data_path 
#    fname_fwd = data_path + meg + 'ico5_forward_5-3L-EMEG-fwd.fif'
#    fname_inv = inv_path + meg + 'InvOp_ico5newreg_fft_1_48_clean_ica_EMEG-inv.fif'
#    vertices_avg = [np.arange(10242), np.arange(10242)]
#    subject_no=subject_inds[cnt]
#
#    forward = mne.read_forward_solution(fname_fwd,surf_ori=True, force_fixed=True)#, surf_ori=True)
#
#    # read inverse operators
#    inverse_operator_eegmeg = read_inverse_operator(fname_inv)
#    
#    snr = 3.0#ediv.mean()
#    lambda2 = 1.0 / snr ** 2
#    method = 'MNE'  # can be 'MNE' or 'sLORETA'
#    mode = 'svd'
#    n_svd_comp = 1
#    lblname_cnt=0
#    sub_src=inverse_operator_eegmeg['src']
#    # compute ctf for all vertices
#    stc_ctf_mne=psf_ctf_28Oct15.cross_talk_function(inverse_operator_eegmeg, forward, labels=None,method=method, lambda2=lambda2, signed=False, mode=mode, n_svd_comp=1, verbose=None)
#    print "psf/ctf finished!"
#
#    subject_from = subjects[subject_no]
#    subject_avg = 'fsaverage'
#    morphmat1=mne.compute_morph_matrix(subject_from=subject_from, subject_to=subject_avg, vertices_from=stc_ctf_mne.vertices, vertices_to=vertices_avg, subjects_dir=data_path)
#    stc_to_pre=stc_ctf_mne.morph_precomputed(subject_avg,vertices_avg,morphmat1, subjects[subject_no])
#    stc_subj_pre = mne.SourceEstimate(stc_to_pre.data[:,:-1].T, vertices=stc_ctf_mne.vertices, tmin=1e-3*2, tstep=1e-3*2, subject=subject_from, verbose=None)
#    stc_to_pre2=stc_subj_pre.morph_precomputed(subject_avg,vertices_avg,morphmat1, subjects[subject_no])
#    stc_to=mne.SourceEstimate(stc_to_pre2.data.T, vertices=vertices_avg, tmin=1e-3*2, tstep=1e-3*2, subject=subject_avg, verbose=None)
#    if cnt==0:
#        stc_all_data=stc_to.data[:,:,np.newaxis]
#    else:
#        stc_all_data=np.concatenate((stc_all_data,stc_to.data[:,:,np.newaxis]),axis=2)
vertices_avg= [np.arange(10242), np.arange(10242)]
subject_avg='fsaverage'
srcin='/imaging/rf02/TypLexMEG/fsaverage/bem/fsaverage_dist-ico-5-src.fif'
#src_avg2=src_avg.copy()
src_avg = mne.read_source_spaces(srcin)
#mne.add_source_space_distances(src_avg)
#stc_mean_data=np.mean(stc_all_data,axis=2)
path_mtx=data_path+'Matx_vert_second.mat'
#scipy.io.savemat(path_mtx,{'stc_mean_data':stc_mean_data})
ltc=scipy.io.loadmat(path_mtx,mat_dtype=True); stc_mean_data=np.squeeze(ltc['stc_mean_data'])
print "stage 1 finished- ctf of split labels"

Matx=stc_mean_data#[:,:-1]# Nv x Nlabels+1
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
## since our numpy version is old have stolen the following bit from newer version for np.unique
assign_vertex[np.logical_or(msort[:,-1]<3, mtokeep<=1)]=-100
aux=assign_vertex.copy()
aux.sort()
flag = np.concatenate(([True], aux[1:] != aux[:-1]))
unver_pre = aux[flag]
idx = np.concatenate(np.nonzero(flag) + ([assign_vertex.size],))
unvercount= np.diff(idx)
#unver=unver_pre[unvercount>=10]
unver_remove=unver_pre[unvercount<=1]
assign_vertex[np.asarray([avi for avi,av in enumerate(assign_vertex) if av in unver_remove])]=-100
#for avii in range(len(assign_vertex)):
#    if assign_vertex[avii] in unver_remove:
#        assign_vertex[avii]=-100
unver=np.unique(assign_vertex)
print len(unver)
#    vcolor=nc(len(unver[:-1]),bytes_=False, cmap='gist_ncar')
labelsv=list()
#num_verts=np.zeros(unver[:-1].shape)
for flc in range(len(unver)-1):
    print flc
    thislabel=np.where(assign_vertex==unver[flc+1])[0]
#    num_verts[flc]=thislabel.shape[0]
    thislabel_data=np.zeros((assign_vertex.shape[0],1))
    thislabel_data[thislabel]=1.0
    thislabel_stc = mne.SourceEstimate(thislabel_data, vertices=vertices_avg, tmin=1e-3*2, tstep=1e-3*2, subject=subject_avg, verbose=None)
    label_append=[x for x in mne.stc_to_label(thislabel_stc,src=src_avg, smooth=True,connected=False, subjects_dir=data_path) if x is not None]
    labelsv.append(label_append)
    if label_append[0].hemi=='lh':#unver[flc+1]<len(assign_vertex)/2.0:
        labelsv[flc][0].name='L'+str(unver[flc+1])+'-lh'
    else:
        labelsv[flc][0].name='L'+str(int(unver[flc+1]-len(assign_vertex)/2))+'-rh'     

    #labelsv[flc][0].color=vcolor[flc]    
labelsvv=[label[0] for label in labelsv]
labels=copy.deepcopy(labelsvv)
num_verts=np.asarray([label.vertices.shape[0] for label in labels])
mean_vert=num_verts.mean()
mean_vert_pre=0.0
mean_vert_change=(mean_vert-mean_vert_pre)/mean_vert
whilec=-1
labnew=list()
labold=[label.name for label in labels]
lendiff=len(labels)

##save onlyverts parcellation
colors=nc(len(labels),bytes_=True, cmap='gist_ncar')
colors=colors[np.random.permutation(len(labels))]
fname='onlyvertices'+'_'+str(len(labels))
make_parcellation_mine_vert(src_avg,labels,colors,subject_avg,fname)
    
while mean_vert_change>0.05:# or lendiff>=5:
    whilec=whilec+1
    print whilec
    lenold=len(labels)
    for cnt, meg in enumerate(ll):
        
        print "subject",cnt
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
        n_svd_comp = 1
        lblname_cnt=0
        sub_src=inverse_operator_eegmeg['src']
#        import time
#        tic=time.time()
        labold=[label.name for label in labels]
        morphed_labels=[label.morph(subject_from='fsaverage', smooth=5, subject_to=subjects[subject_no], subjects_dir=data_path) for label in labels]
#        toc=time.time()
#        print toc-tic
#        morphed_labels=list()
#        for jj,label in enumerate(labels):
#            print jj
#            morphed_labels.append(label.morph(subject_from='fsaverage', smooth=2, grade=vertices_avg, subject_to=subjects[subject_no], subjects_dir=data_path))
        try:
            stc_ctf_mne=psf_ctf_28Oct15.cross_talk_function(inverse_operator_eegmeg, forward, labels=morphed_labels,method=method, lambda2=lambda2, signed=False, mode=mode, n_svd_comp=1, verbose=None)
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
            
                
       
    print "stage 2 finished- loop", whilec
    vertices_avg= [np.arange(10242), np.arange(10242)]
    subject_avg='fsaverage'
    srcin='/imaging/rf02/TypLexMEG/fsaverage/bem/fsaverage_dist-ico-5-src.fif'
    src_avg = mne.read_source_spaces(srcin)
    #mne.add_source_space_distances(src_avg)
    stc_mean_data=np.mean(stc_all_data,axis=2)
    Matx=stc_to.data[:,:-1]# Nv x Nlabels+1
    Matx_zscore=np.zeros(Matx.shape)
    Matx_zscore=zscore(Matx,axis=1)
    assign_vertex=np.argmax(Matx_zscore, axis=1)
    msort=np.sort(Matx_zscore,axis=1)
    mtokeep=msort[:,-1]-msort[:,-2]
    for avii in range(len(assign_vertex)):
        if Matx[avii,assign_vertex[avii]]<np.max(Matx[:,assign_vertex[avii]])/2.:
            assign_vertex[avii]=-100
    
    assign_vertex[np.logical_or(msort[:,-1]<3, mtokeep<=1.5)]=-100
    merge_candidates=np.where(np.logical_and(np.logical_and(msort[:,-1]>=3,msort[:,-2]>=3), mtokeep<=1.5))[0]
    margsort=np.argsort(Matx_zscore,axis=1)
    merge_labels1=margsort[merge_candidates,-2:]
    merge_labels=np.sort(merge_labels1,axis=1)
    new_labels=[tuple(row) for row in merge_labels]
    newlabels=np.unique(new_labels)
    labels_merged=list()
    lblname_cnt=-1
    ## for new merged labels
#    for labcnt in range(newlabels.shape[0]):
##        print labcnt
#        thislabel_verts=np.where(np.logical_and(merge_labels[:,0]==newlabels[labcnt,0],merge_labels[:,1]==newlabels[labcnt,1]))[0]
#        thislabel=merge_candidates[thislabel_verts]
#        if thislabel.shape[0]>=50.0 and labels[newlabels[labcnt,1]].hemi==labels[newlabels[labcnt,0]].hemi:
#            
#            thislabel_data=np.zeros((assign_vertex.shape[0],1))
#            thislabel_data[thislabel]=1.0
#            thislabel_stc = mne.SourceEstimate(thislabel_data, vertices=vertices_avg, tmin=1e-3*2, tstep=1e-3*2, subject=subject_avg, verbose=None)
#            label_append=[x for x in mne.stc_to_label(thislabel_stc,src=src_avg, smooth=True,connected=False, subjects_dir=data_path) if x is not None]
#            if label_append[0].hemi==labels[newlabels[labcnt,1]].hemi:
#                lblname_cnt=lblname_cnt+1
#                labels_merged.append(label_append[0])
#                labels_merged[lblname_cnt].name=labels[newlabels[labcnt,1]].name[:-3]+'_'+labels[newlabels[labcnt,0]].name
    
    ## for the labels we already had
    unver=np.unique(assign_vertex)        
    labvernum=np.zeros(unver.shape[0])
    for uc in range(len(unver)):
        labvernum[uc]=np.where(assign_vertex==unver[uc])[0].shape[0]
        if labvernum[uc]<1.0:
            assign_vertex[np.where(assign_vertex==unver[uc])[0]]=-100
    unver=np.unique(assign_vertex)
    #print unver
    labvernum=np.zeros(unver.shape[0])
    for uc in range(len(unver)):
        labvernum[uc]=np.where(assign_vertex==unver[uc])[0].shape[0]        
    #print sorted(labvernum)
    print labvernum.shape
    lcnt2=-1
    for lcnt in range(len(unver)-1):
#        print lcnt
        thislabel=np.where(assign_vertex==unver[lcnt+1])[0]
        if thislabel.shape[0]>=1.0:
            thislabel_data=np.zeros((assign_vertex.shape[0],1))
            thislabel_data[thislabel]=1.0
            thislabel_stc = mne.SourceEstimate(thislabel_data, vertices=vertices_avg, tmin=1e-3*2, tstep=1e-3*2, subject=subject_avg, verbose=None)
            label_append=[x for x in mne.stc_to_label(thislabel_stc,src=src_avg, smooth=True,connected=False, subjects_dir=data_path) if x is not None]
            if label_append[0].hemi==labels[unver[lcnt+1]].hemi:
                lcnt2=lcnt2+1
                labels_merged.append(label_append[0])
                labels_merged[lblname_cnt+lcnt2+1].name=labels[unver[lcnt+1]].name
    labels=copy.deepcopy(labels_merged)#[label[0] for label in labels_merged]
    num_verts=np.asarray([label.vertices.shape[0] for label in labels])#np.zeros((len(labels),1))
#    for labc in range(len(labels)):
#        num_verts[labc]=labels[labc].vertices.shape[0]
    mean_vert_pre=mean_vert.copy()
    mean_vert=num_verts.mean()
    mean_vert_change=(mean_vert-mean_vert_pre)/mean_vert
    print len(labels)
    lennew=len(labels)
    lendiff=np.abs(lenold-lennew)
    labnew=[label.name for label in labels]
### end of the while loop
    labels_final=[label for label in labels if label.name[-2:]==label.hemi]
    colors=nc(len(labels_final),bytes_=True, cmap='gist_ncar')
    colors=colors[np.random.permutation(len(labels_final))]
    fname='finalpostparc_smooth_stddiff_vertices1w'+str(whilec)+'_'+str(lennew)
    make_parcellation_mine_vert(src_avg,labels_final,colors,subject_avg,fname)
lv2=[len(label.vertices) for label in labels_final]

rs=nc(len(labels_final),bytes_=False, cmap='gist_ncar')
rs=rs[np.random.permutation(len(labels_final))]
for colorc in range(len(labels_final)):
    labels_final[colorc].color=rs[colorc]
mne.write_labels_to_annot(labels_final, subject=subject_avg, parc=subject_avg+'_vertices50n'+str(whilec), subjects_dir=data_path,colormap='gist_ncar',hemi='both',overwrite=True)
lhlabels=[label for label in labels_final if label.hemi=='lh']    
#    labelsff=[label[0] for label in labelsf]
#    mne.write_labels_to_annot(labelsff, subject=subject_from, parc=subject_from+'_finalpostparc_nonsmooth', subjects_dir=data_path,colormap='gist_ncar',hemi='both',overwrite=True)
#    path_c=data_path+'colormatrix.mat'
#    scipy.io.savemat(path_c,{'rs':rs})
#    lfname=[labelsff[ii].name for ii in range(len(labelsff))]
#    lfname_path=data_path+subject_from+'/labels_modified_names.mat'
#    scipy.io.savemat(lfname_path,{'lfname':lfname})
    
#colors=nc(len(labelsvv),bytes_=True, cmap='gist_ncar')
#colors=colors[np.random.permutation(len(labelsvv))]
#fname='onlyvertices10'
#make_parcellation_mine_vert(src_avg,labelsvv,colors,subject_avg,fname)