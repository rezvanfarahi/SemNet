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
#sys.path.insert(1,'/imaging/rf02/TypLexMEG/mne_python0.9')
sys.path.insert(1,'/imaging/local/software/mne_python/latest_v0.9full')

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
from mne.minimum_norm import (read_inverse_operator, make_inverse_operator, write_inverse_operator, apply_inverse, apply_inverse_epochs,psf_ctf_new, point_spread_function)
from mne.minimum_norm.inverse import (prepare_inverse_operator)
import scipy.io
from mne import find_events
from mne.connectivity import seed_target_indices, spectral_connectivity
from mne import read_labels_from_annot

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
    fname_fwd = data_path + meg + 'forward_5-3L-EMEG-fwd.fif'
    fname_inv = inv_path + meg + 'InverseOperator_fft_1_48_clean_ica_EMEG-inv.fif'
    vertices_to = [np.arange(10242), np.arange(10242)]
    subject_no=cnt#subject_inds[0]
#    raw_fname = data_path + meg + 'semloc_ssstf_fft_1_48_clean_ica_raw.fif' #clean_ssp_
#    event_fname = event_path + meg + 'semloc_raw_ssstnew.txt'
#    subjects_dir = data_path 
#
#    tmin = -0.5
#    tmax = 0.7  # Use a lower tmax to reduce multiple comparisons
#
#	#   Setup for reading the raw data
#    raw = io.Raw(raw_fname)
#    events = mne.read_events(event_fname)
#    stim_delay = 0.034 # delay in s
#    event_id={'cncrt_wrd': 1, 'abs_wrd': 2}
#    events[:,0] += np.round( raw.info['sfreq']*stim_delay )   
#	###############################################################################
#	# Read epochs for all channels, removing a bad one
#    
#    picks = mne.pick_types(raw.info, eeg=True, meg=True, eog=False, exclude='bads')
#    if subject_no==3:
#        reject = dict(eeg=150e-6, grad=200e-12, mag=4e-12) #eog=150e-6,
#    else:
#        reject = dict(eeg=120e-6, grad=200e-12, mag=4e-12) #eog=150e-6,
#    epochs = mne.Epochs(raw, events, event_id, tmin, tmax, picks=picks, baseline=(None, 0), proj=True, reject=reject, preload=True)
#    print "epochs finished"
#    epochs.resample(200)
#    epdata=epochs.get_data()
#    em=np.abs(epdata)
#    #emax=np.mean(em[:,:,600:900],axis=2)#
#    sa=em[:,:,110:190]#.max(axis=2)
#    #bm=np.sqrt(np.mean(epdata[:,:,100:500]**2,axis=2))
#    bv=np.var(em[:,:,20:100],axis=2)
#    edv=sa.copy()
#    for ii in range(sa.shape[2]):
#        edv[:,:,ii]=(sa[:,:,ii]**2)/bv
#    
#    ediv=np.mean(np.sqrt(edv),axis=2)
#    print ediv.mean()
    noise_cov = 0#mne.compute_covariance(epochs, method='shrunk', tmin=None, tmax=0.)
    print "noise cov finished!"
    data_cov = 0#mne.compute_covariance(epochs, tmin=0.05, tmax=0.55, method='shrunk')
    print "data cov finished!"
    # read forward solution (sources in surface-based coordinates)
    forward = mne.read_forward_solution(fname_fwd,surf_ori=True, force_fixed=True)#, surf_ori=True)
    forward1 = mne.read_forward_solution(fname_fwd,surf_ori=True, force_fixed=True)#, surf_ori=True)

    # read inverse operators
    inverse_operator_eegmeg = read_inverse_operator(fname_inv)
    #inverse_operator_meg = read_inverse_operator(fname_inv_meg)
    

    labels = [labelss[jj].morph(subject_from='fsaverage', subject_to=subjects[subject_no], subjects_dir=data_path) for jj in range(len(labelss))]

    # regularisation parameter
    snr = 3.0#ediv.mean()
    lambda2 = 1.0 / snr ** 2
    method = 'MNE'  # can be 'MNE' or 'sLORETA'
    mode = 'svd'
    n_svd_comp = 1
    lblname_cnt=0
    whilecnt=0
    while lblname_cnt>-1: 
        whilecnt=whilecnt+1
        print whilecnt
        #mne.io.make_eeg_average_ref_proj(forward['info'])
        #inverse_operator_eegmeg['info']['projs']=inverse_operator_eegmeg['projs']
        #stc_ctf_mne= cross_talk_function(inverse_operator_eegmeg,forward, forward1, noise_cov,data_cov,labels,method=method, mode=mode,signed=False)#, n_svd_comp=n_svd_comp)#
        stc_ctf_mne=psf_ctf_new.cross_talk_function(inverse_operator_eegmeg, forward, labels,method=method, lambda2=lambda2, signed=False, mode=mode, n_svd_comp=1, verbose=None)
        print "psf/ctf finished!"
        #stc_ctf_mne,evoked_fwd = point_spread_function(inverse_operator_eegmeg, forward, labels,method=method, lambda2=lambda2,mode=mode,n_svd_comp=n_svd_comp)#signed=False, 
    
        #stc_psf_meg, _ = point_spread_function(inverse_operator_meg, forward, method=method,labels=labels, lambda2=lambda2, pick_ori='normal', mode=mode,n_svd_comp=n_svd_comp)
    
        # save for viewing in mne_analyze in order of labels in 'labels'
        # last sample is average across PSFs
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
        assign_vertex[logical_or(msort[:,-1]<3, mtokeep<=1)]=-100
        merge_candidates=np.where(logical_and(logical_and(msort[:,-1]>=3,msort[:,-2]>=3), mtokeep<=1))
        margsort=np.argsort(Matx_zscore,axis=1)
        merge_labels1=margsort[merge_candidates[0],-2:]
        merge_labels=np.sort(merge_labels1,axis=1)
        new_labels=[tuple(row) for row in merge_labels]
        newlabels=np.unique(new_labels)
        labels_merged=list()
        lblname_cnt=-1
        for labcnt in range(newlabels.shape[0]):
            print labcnt
            thislabel=np.where(logical_and(merge_labels[:,0]==newlabels[labcnt,0],merge_labels[:,1]==newlabels[labcnt,1]))[0]
            if thislabel.shape[0]>=10:
                lblname_cnt=lblname_cnt+1
                thislabel_data=np.zeros((assign_vertex.shape[0],1))
                thislabel_data[thislabel]=1.0
                thislabel_stc = mne.SourceEstimate(thislabel_data, vertices=stc_ctf_mne.vertices, tmin=1e-3*2, tstep=1e-3*2, subject=subject_from, verbose=None)
                label_append=[x for x in mne.stc_to_label(thislabel_stc,src=subject_from, smooth=False,connected=False, subjects_dir=data_path) if x is not None]
                labels_merged.append(label_append)
                labels_merged[lblname_cnt][0].name=labels[newlabels[labcnt,1]].name[:-3]+'_'+labels[newlabels[labcnt,0]].name
        unver=unique(assign_vertex)
        
        labvernum=np.zeros(unver.shape[0])
        for uc in range(len(unver)):
            labvernum[uc]=np.where(assign_vertex==unver[uc])[0].shape[0]
            if labvernum[uc]<10:
                assign_vertex[np.where(assign_vertex==unver[uc])[0]]=-100
        unver=unique(assign_vertex)
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
            thislabel_stc = mne.SourceEstimate(thislabel_data, vertices=stc_ctf_mne.vertices, tmin=1e-3*2, tstep=1e-3*2, subject=subject_from, verbose=None)
            label_append=[x for x in mne.stc_to_label(thislabel_stc,src=subject_from, smooth=False,connected=False, subjects_dir=data_path) if x is not None]
            labels_merged.append(label_append)
            labels_merged[lblname_cnt+lcnt+1][0].name=labels[unver[lcnt+1]].name
        labels=[label[0] for label in labels_merged]
    
    
    
    # for each individual
#    bb=stc_ctf_mne.in_label(labels)
#    cc=np.abs(stc_ctf_mne.data[:,0])
#    dd=cc.argsort()[::-1]
#    bbv=bb.vertices[0]
#    bbv=np.searchsorted(stc_ctf_mne.vertno[0], bbv)
#    ee=dd[:bb.shape[0]]
#    ff=np.intersect1d(ee,bbv)
    #stc_to.save(data_out)
    # for average
    #stc_to=stc_ctf_mne
    
    ################
#    Matx=stc_to.data
#    stc_all[cnt]=stc_to
#    #Mat_avg=Matx.mean(axis=2)
#    ### combine labels
#    overlap_sub=np.zeros((len(labelss),len(labelss)))
#    for kk, lblk in enumerate(labelss):
#        print kk
#        cc=np.abs(Matx[:,kk])
#        best_vert=sorted(cc,reverse=True)[0]
#        best_verts=np.where(cc>best_vert/2.)[0]
#        
#        for kk2, lblk2 in enumerate(labelss):
#            bb=stc_to.in_label(labelss[kk2])        
#            if kk2%2==0:
#                bbv=bb.vertices[0]
#            else:
#                bbv=bb.vertices[1]+len(stc_to.lh_vertno)
#            
#            overlap_sub[kk,kk2]=len(np.intersect1d(best_verts,bbv))
#    overlap_sub_winner[:,cnt]=np.argmax(overlap_sub, axis=1)
            
    #print rat[:,cnt].mean()
#from scipy.stats import mode
#vertex_winner=mode(assign_vertex, axis=1)
##label_winner=mode(overlap_sub_winner, axis=1)
#path_c=data_path+'aparc_24sep15_vertex_winner.mat'
##path_c2=data_path+'aparc_14sep15_label_winner.mat'
#
#scipy.io.savemat(path_c,{'vertex_winner':vertex_winner})
##scipy.io.savemat(path_c2,{'label_winner':label_winner})
#
#matx_stc = mne.SourceEstimate(vertex_winner[0], vertices=vertices_to,tmin=1e-3 * 1, tstep=1e-3 * 1, subject='fsaverage')
#path3=data_path+'24sep15_ctf_vertex_winner'
#matx_stc.save(path3)
#unver=unique(assign_vertex[:,0])
#labvernum=np.zeros(unver.shape[0])
#for uc in range(len(unver)):
#    labvernum[uc]=np.where(assign_vertex[:,0]==unver[uc])[0].shape[0]
    
#rat_new=rat.copy()
#rat_blnc[::2,:]=np.maximum(rat[::2,:],rat[1::2,:])
#rat_blnc[1::2,:]=np.maximum(rat[::2,:],rat[1::2,:])
#rat_new=rat_blnc.mean(axis=1)
#labl_new=[labl for ii,labl in enumerate(labelss) if rat_new[ii]>=rat_new.mean()-rat_new.std()/2.]#rat_new.mean()]#-rat_new.std()] 
#annot_fname=data_path+'fsaverage/label/rezvan_parc_cnew_mne_ctf_mn.annot'   
#mne.label.write_labels_to_annot(labl_new, subject='fsaverage', parc='rezvan_parc_cnew_mne_ctf_mn', overwrite=True,
#                          subjects_dir=data_path, hemi='both', verbose=None)
#path_rat=data_path+'icaanalysis_results/mne_ratio_cnew2_mat.mat'
#scipy.io.savemat(path_rat,{'rat':rat})
#Matxl=Matxl.transpose()
#Matxr=Matxr.transpose()
#t_Matxl= mne.stats.ttest_1samp_no_p(Matxl, sigma=0, method='relative')
#t_Matxr= mne.stats.ttest_1samp_no_p(Matxr, sigma=0, method='relative')
#t_Matx=np.zeros((20484,2))
#t_Matx[:,0]=t_Matxl
#t_Matx[:,1]=t_Matxr
#fsave_vertices = [np.arange(10242), np.arange(10242)]
#t_stc = mne.SourceEstimate(t_Matx, vertices=fsave_vertices, tmin=1e-3*2, tstep=1e-3*2, subject='fsaverage', verbose=None)
#out_filet=data_path + 'UVttest_CTF_' + lname 
#t_stc.save(out_filet)
#Ml=Matxl.mean(axis=0)
#Mr=Matxr.mean(axis=0)
#M=np.vstack((Ml,Mr)).transpose()
#t_stc = mne.SourceEstimate(M, vertices=fsave_vertices, tmin=1e-3*2, tstep=1e-3*2, subject='fsaverage', verbose=None)
#out_filet=data_path + 'GrandAverage_CTF_' + lname 
#t_stc.save(out_filet)
#	# stc_psf_meg.save('psf_meg')
#"""
#	from mayavi import mlab
#	fmin = 0.
#	time_label = "EEGMEG %d"
#	fmax = stc_psf_eegmeg.data[:, 0].max()
#	fmid = fmax / 2.
#	brain_eegmeg = stc_psf_eegmeg.plot(surface='inflated', hemi='rh',
#		                           subjects_dir=subjects_dir,
#		                           time_label=time_label, fmin=fmin,
#		                           fmid=fmid, fmax=fmax,
#		                           figure=mlab.figure(size=(500, 500)))
#
#	time_label = "MEG %d"
#	fmax = stc_psf_meg.data[:, 0].max()
#	fmid = fmax / 2.
#	brain_meg = stc_psf_meg.plot(surface='inflated', hemi='rh',
#		                     subjects_dir=subjects_dir,
#		                     time_label=time_label, fmin=fmin,
#		                     fmid=fmid, fmax=fmax,
#		                     figure=mlab.figure(size=(500, 500)))
#
#	# The PSF is centred around the right auditory cortex label,
#	# but clearly extends beyond it.
#	# It also contains "sidelobes" or "ghost sources"
#	# in middle/superior temporal lobe.
#	# For the Aud-RH example, MEG and EEGMEG do not seem to differ a lot,
#	# but the addition of EEG still decreases point-spread to distant areas
#	# (e.g. to ATL and IFG).
#	# The chosen labels are quite far apart from each other, so their PSFs
#	# do not overlap (check in mne_analyze)
#"""
