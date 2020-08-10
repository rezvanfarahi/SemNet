"""
=========================================================
TF representation of SemLoc for frequency bands in source labels
=========================================================

"""
# Authors: Alexandre Gramfort <gramfort@nmr.mgh.harvard.edu>
#
# License: BSD (3-clause)


print(__doc__)

import sys
sys.path.append('/imaging/local/software/python_packages/nibabel/1.3.0')
sys.path.append('/imaging/local/software/python_packages/pysurfer/v0.3.1')
sys.path.append('/imaging/local/software/python_packages/scikit-learn/v0.14.1')
sys.path.append('/imaging/local/software/mne_python/latest_v0.11')

# for qsub
# (add ! here if needed) 
sys.path.append('/imaging/local/software/anaconda/latest/x86_64/bin/python')
sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
sys.path.insert(1,'/imaging/local/software/anaconda/2.4.1/3/lib/python3.5/site-packages/scipy')
###

import numpy as np
import mne
import scipy.io as scio

###############################################################################
out_path = '/home/rf02/rezvan/test1/step_by_step/dcm/' # root directory for your MEG data
subjects_dir = '/imaging/rf02/TypLexMEG/'    # where your MRI subdirectories are
# where event-files aremean
data_path = '/imaging/rf02/TypLexMEG/'    # where event files are

label_path = '/home/rf02/rezvan/TypLexMEG/createdlabels_SL/'

inv_path = '/imaging/rf02/TypLexMEG/'

# get indices for subjects to be processed from command line input
# 
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
'meg11_0147/110603/', 


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
labelss = mne.read_labels_from_annot(subject='fsaverage', parc='rez_aparc_24sep15', subjects_dir=data_path)
#labelss.remove(labelss[-3])
#labelss.remove(labelss[-2])
labelss.remove(labelss[-1])
for lbl in labelss:
    lbl.values.fill(1.0)
label_names=[label.name for label in labelss]
label_index=np.array([label_names.index('ATL_latmed-lh'),label_names.index('supramarginal-lh'),label_names.index('posterior_inferiorfrontal-lh'),label_names.index('middle_temporal-lh')])#label_names.index('ATL_latmed-lh'),

labellist=[labelss[ii] for ii in label_index]
#labellist = ['lh.ATL_latmed.label']#['atlleft-lh.label']#['lh.ATL_latmed.label']
labellist2 = ['lh.ATL_DCM-lh.label','lh.WFA_DCM-lh.label']
#lfname=data_path+'fsaverage/label/'+labellist2[0]#label_path+labellist[0]#data_path+'fsaverage/label/'+labellist[0]
#label0=mne.read_label(lfname,subject='fsaverage')
#labellist.append(label0)
lfname=data_path+'fsaverage/label/'+labellist2[1]#label_path+labellist[0]#data_path+'fsaverage/label/'+labellist[0]
label1=mne.read_label(lfname,subject='fsaverage')
labellist.append(label1)
print labellist
#labellist = ['lh.WFA_DCM-lh.label','lh.ATL_DCM-lh.label']

# Labels/ROIs
#labellist = ['atlleft-lh']#, 'atlleftt-rh'
tmin, tmax = -0.5, 0.7
# frequency bands with variable number of cycles for wavelets
stim_delay = 0.034 # delay in s

tmin1=-500
tstep1=700
stc_allc=range(17)
vertices_to = [np.arange(10242), np.arange(10242)]
stc_alla=range(17)
vertices_to = [np.arange(10242), np.arange(10242)]
cl=range(5)
al=range(5)
clw=np.arange(5)
alw=np.arange(5)
mlw=np.arange(5)
cncrt_mat=np.zeros((5,1201))
cons=np.zeros((17,5,5))
abs_mat=np.zeros((5,1201))
for thislabelc,thislabel in enumerate(labellist):
    for ii, meg in enumerate(ll):
        print ii
        srcin='/imaging/rf02/TypLexMEG/fsaverage/bem/fsaverage-ico-5-src.fif'
        src = mne.read_source_spaces(srcin)
        fname = data_path + meg + 'firstMorphed_ico_signed_typlex_ica_word_Source_Evoked_m500_700' 
        stc_allc = mne.read_source_estimate(fname)
        
        
        cl1=stc_allc.in_label(thislabel).data
    #    cl1w=np.argmax(np.mean(np.abs(cl1[:,550:950]),axis=1))
    #    print cl1w
        
        
        fname = data_path + meg + 'firstMorphed_ico_signed_typlex_ica_nonword_Source_Evoked_m500_700'  
        stc_alla = mne.read_source_estimate(fname)
        al1=stc_alla.in_label(thislabel).data
    #    al1w=np.argmax(np.mean(np.abs(al1[:,550:950]),axis=1))
    #    print al1w
        clm=np.mean((cl1[:,550:750]),axis=1)
        alm=np.mean((al1[:,550:750]),axis=1)
        if ii==0:
            slw=(clm-alm)[np.newaxis,:]
        else:
            slw=np.append(slw,(clm-alm)[np.newaxis,:],axis=0)
    MT,pT,tT= mne.stats.permutation_t_test(slw)#MT= mne.stats.ttest_1samp_no_p(slw, sigma=1e-3, method='relative')
    mlw[thislabelc]=np.argmax(MT)

for ii, meg in enumerate(ll):
    print ii
    srcin='/imaging/rf02/TypLexMEG/fsaverage/bem/fsaverage-ico-5-src.fif'
    src = mne.read_source_spaces(srcin)
    fname = data_path + meg + 'firstMorphed_ico_signed_typlex_ica_word_Source_Evoked_m500_700' 
    stc_allc = mne.read_source_estimate(fname)
#    lfname=data_path+'fsaverage/label/'+labellist[0]
#    label1=mne.read_label(lfname,subject='fsaverage')
#    lfname=data_path+'fsaverage/label/'+labellist[1]
#    label2=mne.read_label(lfname,subject='fsaverage')
    for cntcl in range(len(cl)):        
        cl[cntcl]=stc_allc.in_label(labellist[cntcl]).data
#        clw[cntcl]=np.argmax(np.mean(np.abs(cl[cntcl][:,500:1000]),axis=1))
#    print clw
    fname = data_path + meg + 'firstMorphed_ico_signed_typlex_ica_nonword_Source_Evoked_m500_700'  
    stc_alla = mne.read_source_estimate(fname)
    for cntal in range(len(al)):        
        al[cntal]=stc_alla.in_label(labellist[cntal]).data
#        alw[cntal]=np.argmax(np.mean(np.abs(al[cntal][:,500:1000]),axis=1))
#    print alw
    for cntw in range(len(mlw)):
#        if alw[cntw]!=clw[cntw]:
#            ml1=np.argmax(np.array([np.max(np.mean(np.abs(cl[cntw][:,500:1000]),axis=1)),np.max(np.mean(np.abs(al[cntw][:,500:1000]),axis=1))]))
#            bl1=np.array([clw[cntw],alw[cntw]])
#            clw[cntw]=bl1[ml1]; alw[cntw]=bl1[ml1]
        thislabelsign=mne.label.label_sign_flip(labellist[cntw],src)[mlw[cntw]]
        cncrt_mat[cntw,:]=thislabelsign*np.expand_dims(cl[cntw][mlw[cntw],:].transpose(),0)
        abs_mat[cntw,:]=thislabelsign*np.expand_dims(al[cntw][mlw[cntw],:].transpose(),0)

#    if al2w!=cl2w:
#        ml2=np.argmax(np.array([np.max(np.mean(np.abs(cl[cntw][:,500:1000]),axis=1)),np.max(np.mean(np.abs(al[cntw][:,500:1000]),axis=1))]))
#        bl2=np.array([cl2w,al2w])
#        cl2w=bl2[ml2]; al2w=bl2[ml2]
#    
#    label1sign=mne.label.label_sign_flip(label1,src)[al1w]
#    label2sign=mne.label.label_sign_flip(label2,src)[al2w]
#    cncrt_label1=label1sign*np.expand_dims(cl1[cl1w,:].transpose(),0)#mne.extract_label_time_course(stc_allc, label1, src, mode='mean_flip')
#    cncrt_label2=label2sign*np.expand_dims(cl2[cl2w,:].transpose(),0)#mne.extract_label_time_course(stc_allc, label2, src, mode='mean_flip')
#    cncrt_mat=np.concatenate((cncrt_label1,cncrt_label2),0)
#    abs_label1=label1sign*np.expand_dims(al1[al1w,:].transpose(),0)#mne.extract_label_time_course(stc_alla, label1, src, mode='mean_flip')
#    abs_label2=label2sign*np.expand_dims(al2[al2w,:].transpose(),0)#mne.extract_label_time_course(stc_alla, label2, src, mode='mean_flip')
#    abs_mat=np.concatenate((abs_label1,abs_label2),0)       
    conabs_mat=np.concatenate((np.expand_dims(cncrt_mat,2),np.expand_dims(abs_mat,2)),2)
    con=mne.connectivity.spectral_connectivity(data=cncrt_mat[np.newaxis,:,550:750],method='imcoh',sfreq=1000.,mode='multitaper',fmin=1.,fmax=45.,mt_adaptive=True,faverage=True)
    conc=np.corrcoef(cncrt_mat[:,550:750])#con[0]+con[0].transpose(1,0,2)
    con2=mne.connectivity.spectral_connectivity(data=abs_mat[np.newaxis,:,550:750],method='imcoh',sfreq=1000.,mode='multitaper',fmin=1.,fmax=45.,mt_adaptive=True,faverage=True)
    cona=np.corrcoef(abs_mat[:,550:750])#con2[0]+con2[0].transpose(1,0,2)
    cons[ii,:,:]=np.squeeze(conc-cona)
    
    out_file=out_path + meg[:10] + '_SemLoc_Evoked_5ROIs_forDCM.mat'
#    scio.savemat(out_file,{'conabs_mat':conabs_mat})
consr=np.reshape(cons,(17,25))
import scipy
MTc= mne.stats.ttest_1samp_no_p(consr, sigma=1e-3, method='relative')#Tc,pTc,tTc= mne.stats.permutation_t_test(consr)
pTc=scipy.stats.t.sf(np.abs(MTc),16)*2
pfinal=np.reshape(pTc,(5,5))
print pfinal
