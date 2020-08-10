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
from mne.fixes import in1d, digitize
from mne import Label
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
#    fname_inv = inv_path + meg + 'InvOp_ico5diagreg_fft_1_48_clean_ica_EMEG-inv.fif'
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
assign_vertex[np.logical_or(msort[:,-1]<3, mtokeep<1.0)]=-100
aux=assign_vertex.copy()
aux.sort()
flag = np.concatenate(([True], aux[1:] != aux[:-1]))
unver_pre = aux[flag]
idx = np.concatenate(np.nonzero(flag) + ([assign_vertex.size],))
unvercount= np.diff(idx)
#unver=unver_pre[unvercount>=10]
unver_remove=unver_pre[unvercount<=1]
if unver_remove.any():
    assign_vertex[np.asarray([avi for avi,av in enumerate(assign_vertex) if av in unver_remove])]=-100
#for avii in range(len(assign_vertex)):
#    if assign_vertex[avii] in unver_remove:
#        assign_vertex[avii]=-100
unver=np.unique(assign_vertex)[1:]
unver=unver[unver<10242]
Mdiag=np.diag(Matx_zscore)
Mdiag_unver=Mdiag[unver]
assign_vertex_new=-100*np.ones(assign_vertex.shape)
for rgcnt in range(len(Mdiag_unver)):
    print rgcnt
    if np.min(Mdiag_unver)<1000:
        this_seed=unver[np.argmin(Mdiag_unver)]
        this_seed_realm=np.where(Matx[:,this_seed]>=np.max(Matx[:,this_seed])/2.)[0]
        assign_vertex_new[this_seed_realm]=int(this_seed)
        Mdiag_unver[np.argmin(Mdiag_unver)]=1000
#    ##remove overlaps from unver
#    inur=np.intersect1d(unver,this_seed_realm)
#    for incnt in inur:
#        unver[np.where(unver==incnt)[0]]=-100
#    unver=np.unique(unver)


#path_mtx=data_path+'lefthem_winner_vertices.mat'
#"""
#scipy.io.savemat(path_mtx,{'unver':unver})
#"""
#path_mtx=data_path+'bothhem_winner_vertices.mat'
##scipy.io.savemat(path_mtx,{'stc_mean_data':stc_mean_data})
#ltc=scipy.io.loadmat(path_mtx,mat_dtype=True); unver=np.squeeze(ltc['unver'])
aux=assign_vertex_new.copy()
aux.sort()
flag = np.concatenate(([True], aux[1:] != aux[:-1]))
unver_pre = aux[flag]
idx = np.concatenate(np.nonzero(flag) + ([assign_vertex_new.size],))
unvercount= np.diff(idx)
#unver=unver_pre[unvercount>=10]
unver_remove=unver_pre[unvercount<=25]
if unver_remove.any():
    assign_vertex_new[np.asarray([avi for avi,av in enumerate(assign_vertex_new) if av in unver_remove])]=-100
unver=np.unique(assign_vertex_new)[1:]
print len(unver)
#    vcolor=nc(len(unver[:-1]),bytes_=False, cmap='gist_ncar')
labelsv=list()
#num_verts=np.zeros(unver[:-1].shape)
leftright_labels=range(2*len(unver))
path_mtx='/imaging/rf02/parcellation/MNIlh.mat'
MNIlh_pre=scipy.io.loadmat(path_mtx,mat_dtype=True); MNIlh=np.squeeze(MNIlh_pre['MNIlh'])
path_mtx='/imaging/rf02/parcellation/MNIrh.mat'
MNIrh_pre=scipy.io.loadmat(path_mtx,mat_dtype=True); MNIrh=np.squeeze(MNIrh_pre['MNIrh'])

for flc in range(len(unver)):
    print flc
    thislabel=unver[flc]#np.where(assign_vertex==unver[flc+1])[0]
#    num_verts[flc]=thislabel.shape[0]
    thislabel_left=np.zeros((assign_vertex_new.shape[0],1))
    thislabel_left=np.where(assign_vertex_new==thislabel)[0]
    thislabel_left=thislabel_left[thislabel_left<10242]
    leftright_labels[flc*2]=thislabel_left
    thislabel_right=thislabel_left.copy()
    for thislabii in range(len(thislabel_left)):
        this_vertex=MNIlh[thislabel_left[thislabii],:]
        this_dist=np.sqrt(np.sum((MNIrh-this_vertex)**2,1))
        min_dist=np.min(this_dist)
        w_min_dist=np.argmin(this_dist)
        while w_min_dist+10242 in thislabel_right:#np.sum(thislabel_left==w_min_dist+10242)==1:
            this_dist[w_min_dist]=20000;
            min_dist=np.min(this_dist)
            w_min_dist=np.argmin(this_dist)
        thislabel_right[thislabii]=w_min_dist+10242;
    leftright_labels[flc*2+1]=thislabel_right
    subject_avg='fsaverage'
    srcin='/imaging/rf02/TypLexMEG/fsaverage/bem/fsaverage_dist-ico-5-src.fif'
    src_avg = mne.read_source_spaces(srcin)
for flc in range(len(leftright_labels)):
    print flc
    #    thislabel=unver[flc]#np.where(assign_vertex==unver[flc+1])[0]
    #    num_verts[flc]=thislabel.shape[0]
    thislabel_data=np.zeros((assign_vertex_new.shape[0],1))
    thislabel_data[leftright_labels[flc]]=1.0
    thislabel_stc = mne.SourceEstimate(thislabel_data, vertices=vertices_avg, tmin=1e-3*2, tstep=1e-3*2, subject=subject_avg, verbose=None)
    label_append=[x for x in mne.stc_to_label(thislabel_stc,src=src_avg, smooth=True,connected=False, subjects_dir=data_path) if x is not None]
    label_append[0].smooth(smooth=2,subject='fsaverage',subjects_dir=data_path)
    if label_append[0].hemi=='lh':#unver[flc+1]<len(assign_vertex)/2.0:
        label_append[0].name='L'+str(int(unver[flc/2]))+'-lh'
        #hemi_src = src_avg[0]
    else:
        label_append[0].name='L'+str(int(unver[flc/2]))+'-rh' 
        #hemi_src = src_avg[1]
    labelsv.append(label_append[0])
#    nearest = hemi_src['nearest']
#    # find new vertices
#    include = mne.fixes.in1d(nearest, label_append[0].vertices, False)
#    vertices = np.nonzero(include)[0]
#    # values
#    nearest_in_label = digitize(nearest[vertices], label_append[0].vertices, True)
#    values = label_append[0].values[nearest_in_label]
#    # pos
#    pos = hemi_src['rr'][vertices]
#    label_appendf = Label(vertices, pos, values, label_append[0].hemi, label_append[0].comment, label_append[0].name,
#                      None, label_append[0].subject, label_append[0].color)
#    labelsv.append(label_appendf)
    

        
        
    #labelsv[flc][0].color=vcolor[flc]    
labelsvv=[label for label in labelsv]
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
labels_original= copy.deepcopy(labelsvv)
colors_orig=nc(len(labels_original)/2,bytes_=True, cmap='hsv')
colors_orig=colors_orig[np.random.permutation(len(labels_original)/2)]
#colors=nc(len(labels),bytes_=True, cmap='gist_ncar')
#colors=colors[np.random.permutation(len(labels))]
fname='onlyvertices'+'_actualRG_'+str(len(labelsvv))
make_parcellation_mine(src_avg,labels_original,colors_orig,subject_avg,fname)
labels_original_name = [label.name[:-3] for label in labels if label.name.endswith('lh')]#
#mne.read_labels_from_annot(subject='fsaverage', parc='fsaverage_'+fname, subjects_dir=subjects_dir)#parc='rez_aparc_24sep15'fsaverage_morphed_split_aparc

print "finish now if you wish"
    
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
    
    assign_vertex[np.logical_or(msort[:,-1]<3, mtokeep<1)]=-100
    merge_candidates=np.where(np.logical_and(np.logical_and(msort[:,-1]>=3,msort[:,-2]>=3), mtokeep<=1.0))[0]
    margsort=np.argsort(Matx_zscore,axis=1)
    merge_labels1=margsort[merge_candidates,-2:]
    merge_labels=np.sort(merge_labels1,axis=1)
    new_labels=[tuple(row) for row in merge_labels]
    newlabels=[]
    for new_label in new_labels:
        if new_label not in newlabels:
            newlabels.append(new_label)
    newlabels=np.asarray(newlabels)
    for nlcnt, nl in enumerate(newlabels):
        this_nl=[nl0 for nl0 in new_labels if np.all(np.asarray(nl0)==newlabels[nlcnt])]
        contr_label1=np.where(assign_vertex==newlabels[nlcnt][0])[0]
        contr_label2=np.where(assign_vertex==newlabels[nlcnt][1])[0]
    labels_merged=list()
    lblname_cnt=-1

    ## for new merged labels
#    for labcnt in range(newlabels.shape[0]):
#        #        print labcnt
#        merlabel_verts=np.where(np.logical_and(merge_labels[:,0]==newlabels[labcnt,0],merge_labels[:,1]==newlabels[labcnt,1]))[0]
#        len_merlabel=len(merlabel_verts)
#        contr_label0=np.where(assign_vertex==newlabels[labcnt][0])[0]
#        contr_label1=np.where(assign_vertex==newlabels[labcnt][1])[0]
#        thislabel=merge_candidates[thislabel_verts]
#        if thislabel.shape[0]>=25.0 and labels[newlabels[labcnt,1]].hemi==labels[newlabels[labcnt,0]].hemi:
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
    
    print len(labels)
    lennew=len(labels)
    lendiff=np.abs(lenold-lennew)
    labnew=[label.name for label in labels]
    ### end of the while loop
    labels_final=[label for label in labels if label.name[-2:]==label.hemi]    
    lnames=[label.name for label in labels_final]
    balanced_labels=list()
    for lcnt01 in range(len(labels_final)):
        if lnames[lcnt01].endswith('lh'):
            for lcnt02 in range(len(labels_final)):
                if lnames[lcnt02].endswith('rh') and lnames[lcnt02][:-3]==lnames[lcnt01][:-3] and labels_final[lcnt02].hemi!=labels_final[lcnt01].hemi:
                    balanced_labels.append(labels_final[lcnt01])
                    balanced_labels.append(labels_final[lcnt02])
    print "save the results! ^v^"
    colors_fin=nc(len(balanced_labels)/2,bytes_=True, cmap='gist_ncar')
    colors_fin=colors_fin[np.random.permutation(len(balanced_labels)/2)]
    for labcntc,label in enumerate(balanced_labels[::2]):
        if label.name[:-3] in labels_original_name:
            colors_fin[labcntc,:]=colors_orig[labels_original_name.index(label.name[:-3]),:]
    fname='finalpostparc_smooth_actualRG_leftonly_'+str(whilec)+'_'+str(len(balanced_labels))
    make_parcellation_mine(src_avg,balanced_labels,colors_fin,subject_avg,fname)
    lv2=[len(label.vertices) for label in balanced_labels]
#    blname=[label.name for label in balanced_labels]   
#    colors=nc(len(balanced_labels)/2,bytes_=True, cmap='gist_ncar')
#    colors=colors[np.random.permutation(len(balanced_labels)/2)]
#    fname='finalpostparc_smooth_stddiff_LRmirror_bal_wmerge_'+str(whilec)+'_'+str(len(balanced_labels))
#    make_parcellation_mine(src_avg,balanced_labels,colors,subject_avg,fname)
    labels=copy.deepcopy(balanced_labels)
    num_verts=np.asarray([label.vertices.shape[0] for label in labels])#np.zeros((len(labels),1))
#    for labc in range(len(labels)):
#        num_verts[labc]=labels[labc].vertices.shape[0]
    mean_vert_pre=mean_vert.copy()
    mean_vert=num_verts.mean()
    mean_vert_change=0#(mean_vert-mean_vert_pre)/mean_vert
#rs=nc(len(labels_final),bytes_=False, cmap='gist_ncar')
#rs=rs[np.random.permutation(len(labels_final))]
#labels_fc = mne.read_labels_from_annot(subject='fsaverage', parc='fsaverage_'+fname, subjects_dir=subjects_dir)#parc='rez_aparc_24sep15'fsaverage_morphed_split_aparc
#for label in labels_fc:
#    if label in labels_original:
#        label.color=labels_original[labels_original.index(label)].color
#mne.write_labels_to_annot(labels_fc, subject=subject_avg, parc='fsaverage_scolored_'+fname, subjects_dir=data_path,colormap='gist_ncar',hemi='both',overwrite=True)




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