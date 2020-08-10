# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 17:03:32 2016

@author: rf02
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 14:55:11 2015

@author: rf02
"""

# Author: Martin Luessi <mluessi@nmr.mgh.harvard.edu>
#
# License: BSD (3-clause)

print(__doc__)
import sys

#sys.path.append('/imaging/local/software/python_packages/nibabel/1.3.0')
#sys.path.append('/imaging/local/software/python_packages/pysurfer/v0.4')
#sys.path.append('/imaging/local/software/python_packages/scikit-learn/v0.14.1')
#sys.path.append('/imaging/local/software/mne_python/latest')

# for qsub
# (add ! here if needed)
sys.path.append('/imaging/local/software/EPD/latest/x86_64/bin/python')
sys.path.insert(1,'/imaging/rf02/TypLexMEG/mne_python0.9')
#sys.path.insert(1,'/imaging/rf02/mne_python_11')

#sys.path.insert(1,'/imaging/local/software/mne_python/latest_v0.9')
#sys.path.append('/imaging/local/software/python_packages/scikit-learn/v0.16.1')
sys.path.insert(1,'/imaging/local/software/python_packages/scikit-learn/v0.13.1/lib.linux-x86_64-2.7/sklearn/externals')
sys.path.append('/home/rf02/rezvan/test1')
#sys.path.append('/imaging/local/software/python_packages/scikit-learn/v0.14.1')
#sys.path.insert(1,'/imaging/local/software/python_packages/scikit-learn/v0.16.1')
sys.path.insert(1,'/imaging/rf02/scikit-learn-0.16.1')
sys.path.insert(1,'/home/rf02/.local/lib/python2.7/site-packages')

import joblib
###


import os
import numpy as np
import mne
import scipy
from scipy.stats.mstats import zscore
#import sklearn
from mne.io import Raw
from mne.minimum_norm import (read_inverse_operator, psf_ctf_28Oct15, apply_inverse, apply_inverse_epochs)
from mne.connectivity import seed_target_indices, spectral_connectivity
from numpy import *
import numpy
import time

def entropy(counts):
    '''Compute entropy.'''
    ps = counts/float(sum(counts))  # coerce to float and normalize
    ps = ps[nonzero(ps)]            # toss out zeros
    H = -sum(ps * numpy.log2(ps))   # compute entropy
    
    return H

def mi(x, y, bins):
    '''Compute mutual information'''
    counts_xy = histogram2d(x, y, bins=bins)[0]
    counts_x  = histogram(  x,    bins=bins)[0]
    counts_y  = histogram(  y,    bins=bins)[0]
    
    H_xy = entropy(counts_xy)
    H_x  = entropy(counts_x)
    H_y  = entropy(counts_y)
    
    return H_x + H_y - H_xy
data_path = '/imaging/rf02/TypLexMEG/' # root directory for your MEG data
os.chdir(data_path)
event_path = '/group/erp/data/olaf.hauk/MEG/TypLexM/data/MF22_new/'   # where event files are
label_path = '/imaging/rf02/TypLexMEG/fsaverage/label/'
inv_path = '/imaging/rf02/TypLexMEG/'
subjects_dir = data_path
print sys.argv
subject_inds = [0]
for ss in sys.argv[1:]:
   subject_inds.append( int( ss ) )

#subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
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




#label_name_lh = ['atlleft-lh', 'atlright-rh']


# Labels/ROIs'lh.lateralorbitofrontal','rh.lateralorbitofrontal']#
#labellist = ['lh.ATL_latmed','lh.supramarginal','lh.posterior_inferiorfrontal','lh.middle_temporal']#[ 'lh.posterior_fusiform','rh.posterior_fusiform', 'lh.lateraloccipital','rh.lateraloccipital','lh.precentral','rh.precentral', 'lh.ATL_latmed','rh.ATL_latmed']#labellist = [ 'leftATL-lh','rightATL-rh']#, 'medialtempright-rh','medialtempleft-lh']'lh.postcentral','rh.postcentral',
labelss = mne.read_labels_from_annot(subject='fsaverage', parc='rez_aparc_25oct16', subjects_dir=data_path)
#labelss.remove(labelss[-3])
#labelss.remove(labelss[-2])
#labelss.remove(labelss[-1])
for lbl in labelss:
    lbl.values.fill(1.0)
label_names=[label.name for label in labelss]
label_index=np.array([label_names.index('ATL_latmed-lh'),label_names.index('ATL_latmed-rh'),\
label_names.index('AG_new-lh'),label_names.index('AG_new-rh'),\
label_names.index('posterior_inferiorfrontal-lh'),label_names.index('posterior_inferiorfrontal-rh'),\
label_names.index('middle_temporal-lh'),label_names.index('middle_temporal-rh'),\
label_names.index('lateralorbitofrontal-lh'),label_names.index('lateralorbitofrontal-rh'),\
label_names.index('inferior2_precentral-lh'),label_names.index('inferior2_precentral-rh'),\
label_names.index('superior2_precentral-lh'),label_names.index('superior2_precentral-rh'),\
label_names.index('inferior2_postcentral-lh'),label_names.index('inferior2_postcentral-rh'),\
label_names.index('superior2_postcentral-lh'),label_names.index('superior2_postcentral-rh'),\
label_names.index('lateraloccipital-lh'),label_names.index('lateraloccipital-rh'),\
#label_names.index('vWFA2-lh'),label_names.index('vWFA2-rh'),\
])
labellist1=[labelss[ii] for ii in label_index]
#labellist=labellist1[:10]+[labellist1[10]+labellist1[12]+labellist1[14]+labellist1[16]]+[labellist1[11]+labellist1[13]+labellist1[15]+labellist1[17]]+[labellist1[18]+labellist1[20]]+[labellist1[19]+labellist1[21]]
labellist=labellist1[:10]+[labellist1[10]+labellist1[12]]+[labellist1[11]+labellist1[13]]+[labellist1[14]+labellist1[16]]+[labellist1[15]+labellist1[17]]+labellist1[18:]#[labellist1[18]+labellist1[20]]+[labellist1[19]+labellist1[21]]

print labellist
#label_names.index('medialorbitofrontal-lh'),label_names.index('medialorbitofrontal-rh'),\
#
#mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/lh.supramarginal.label', supramarginal)
#
#posterior_inferiorfrontal=labelss[label_name.index('posterior_inferiorfrontal-lh')]
#mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/lh.posterior_inferiorfrontal.label', posterior_inferiorfrontal)
#
#middle_temporal=labelss[label_name.index('middle_temporal-lh')]
#mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/lh.middle_temporal.label', middle_temporal)



print "stage 1 finished- ctf of split labels"
vertices_avg= [np.arange(10242), np.arange(10242)]
subject_avg='fsaverage'
srcin='/imaging/rf02/TypLexMEG/fsaverage/bem/fsaverage_dist-ico-5-src.fif'
#src_avg2=src_avg.copy()
src_avg = mne.read_source_spaces(srcin)
#mne.add_source_space_distances(src_avg)

    
    
for meg in ll:
    
    srcin='/imaging/rf02/TypLexMEG/fsaverage/bem/fsaverage-ico-5-src.fif'
    src = mne.read_source_spaces(srcin)#inverse_operator['src']  # the source space used

    stc_incoh = data_path + meg + 'firstMorphed_ico_SemLoc_ica_Abstract_Source_Evoked_m500_700'
    stc_coh = mne.read_source_estimate(stc_incoh)
    #stc_coh.resample(100)
    
    seed_ts = mne.extract_label_time_course(stc_coh, labellist, src, mode='pca_flip')
    
# Restrict the source estimate to the label in the left auditory cortex
    #stcp = apply_inverse(evoked, inverse_operator, lambda2, method, pick_ori="normal")
    #stc = mne.morph_data(subject_from, subject_to, stcp, subjects_dir=data_path, n_jobs=1, grade=vertices_to)
#        stc_label = stc_coh.in_label(label1)
#
## Find number and index of vertex with most power
#        src_pow = np.sum(stc_label.data ** 2, axis=1)
#        if label_name[0] is 'l':
#            print label_name
#            seed_vertno = stc_label.vertices[0][np.argmax(src_pow)]
#            seed_idx = np.searchsorted(stc_coh.vertices[0], seed_vertno)  # index in original stc
#            print seed_idx
#
#        else:
#            print label_name
#            seed_vertno = stc_label.vertices[1][np.argmax(src_pow)]
#            seed_idx = stc_coh.vertices[0].shape[0]+ np.searchsorted(stc_coh.vertices[1], seed_vertno)  # index in original stc
#            print seed_idx

# Generate index parameter for seed-based connectivity analysis
#        n_sources = stc_coh.data.shape[0]
#        indices = seed_target_indices([0], np.arange(n_sources))
#    vertices = [src[i]['vertno'] for i in range(2)]
#    n_signals_tot = 1 + len(vertices[0]) + len(vertices[1])
#
#    indices = seed_target_indices([0], np.arange(1,n_signals_tot))
#
#
## Compute inverse solution and for each epoch. By using "return_generator=True"
## stcs will be a generator object instead of a list. This allows us so to
## compute the Coherence without having to keep all source estimates in memory.
#
## Now we are ready to compute the Coherence in the alpha and beta band.
## fmin and fmax specify the lower and upper freq. for each band, resp.
#    fmin = (4., 8., 13., 31.)
#    fmax = (7., 12., 30, 45.)
#    tmin2 = 0.0
#    tmax2 = 0.551
#    tminb2 = -0.4
#    tmaxb2 = 0.
#    sfreq = raw.info['sfreq']  # the sampling frequency
#    #cwt_n_cycles=np.array([2,2,2,2,3,3,3,3,3,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7])
#    cwt_frequencies = np.arange(4, 46, 1)#np.array([6,10,22,38])#
#    cwt_n_cycles = cwt_frequencies / float(3)#np.array([2,3,5,7])#
#    #cwt_n_cycles[cwt_frequencies<=20] = 2
#    #cwt_n_cycles[cwt_frequencies<=20] += np.arange(0,1,0.059)#0.117)#0.059)
#    #con_blcn, freqs_blcn, times_blcn, n_epochs_blcn, n_tapers_blcn = spectral_connectivity(label_ts_Abstract, method=con_methods, mode='multitaper', sfreq=sfreq, fmin=fmin, fmax=fmax, faverage=True, tmin=tminb2, tmax=tmaxb2, mt_adaptive=True, n_jobs=2)
#
#    #con_blab, freqs_bl, times_bl, n_epochs_bl, n_tapers_bl = spectral_connectivity(label_ts_abstract, method=con_methods, mode='multitaper', sfreq=sfreq, fmin=fmin, fmax=fmax, faverage=True, tmin=tminb2, tmax=tmaxb2, mt_adaptive=True, n_jobs=2)
#
#
## Now we compute connectivity. To speed things up, we use 2 parallel jobs
## and use mode='multitaper', mt_adaptive=True, which uses a FFT with a Hanning window
## to compute the spectra (instead of multitaper estimation, which has a
## lower variance but is slower). By using faverage=True, we directly
## average the Coherence in the alpha and beta band, i.e., we will only
## get 2 frequency bins
#	
#    print "spectral connectivity coherence"
#    coh=np.zeros((len(labellist),len(labellist),4,2))
#    for actc,tc in enumerate(range(150,251,100)): #time count
#        tminw=tc/1000.
#        tmaxw=(tc+200)/1000.
#        coh[:,:,:,actc], freqs1, times1, n_epochs1, _ = spectral_connectivity(seed_ts, method='coh', mode='multitaper', sfreq=sfreq, fmin=fmin, fmax=fmax, faverage=True, tmin=tminw, tmax=tmaxw, n_jobs=4)
###spectral_connectivity(data, method='coh', indices=None, sfreq=6.283185307179586, mode='multitaper', fmin=None, fmax=inf, fskip=0, faverage=False, tmin=None, tmax=None, mt_bandwidth=None, mt_adaptive=False, mt_low_bias=True, cwt_frequencies=None, cwt_n_cycles=7, block_size=1000, n_jobs=1, verbose=None)
#    print('coh fin ')
#    import scipy.io as scio
#    path_mtx=data_path + meg +'coh_14ROIs_2_abstract.mat'
#    scio.savemat(path_mtx,{'coh':coh})
#
## Generate a SourceEstimate with the Coherence. This is simple since we
## used a single seed. For more than one seeds we would have to split coh.
## Note: We use a hack to save the frequency axis as time
#    print "use a hack to save the frequency axis as time"
##    tmin1 = 0.150#np.mean(freqs[0])
##    tstep1 = 0.100#np.mean(freqs[1]) - tmin
##    coh_stc = mne.SourceEstimate(coh[:,0,:], vertices=stc_coh.vertices, tmin= tmin1, tstep= tstep1, subject='fsaverage')
##    data_out = data_path + meg +  'Morphed_ico_SemLoc_ica_equalized_Abstract_Mlttap_new_ctf_Coherence_Theta_150_450_200ms_' + thislabelname
##    print 'theta'
##    coh_stc.save(data_out)
##    
##    coh_stc = mne.SourceEstimate(coh[:,1,:], vertices=stc_coh.vertices, tmin= tmin1, tstep= tstep1, subject='fsaverage')
##    data_out = data_path + meg +  'Morphed_ico_SemLoc_ica_equalized_Abstract_Mlttap_new_ctf_Coherence_Alpha_150_450_200ms_' + thislabelname
##    print 'alpha'
##    coh_stc.save(data_out)
##    
##    coh_stc = mne.SourceEstimate(coh[:,2,:], vertices=stc_coh.vertices, tmin= tmin1, tstep= tstep1, subject='fsaverage')
##    data_out = data_path + meg +  'Morphed_ico_SemLoc_ica_equalized_Abstract_Mlttap_new_ctf_Coherence_Beta_150_450_200ms_' + thislabelname
##    print 'beta'
##    coh_stc.save(data_out)
##    
##    coh_stc = mne.SourceEstimate(coh[:,3,:], vertices=stc_coh.vertices, tmin= tmin1, tstep= tstep1, subject='fsaverage')
##    data_out = data_path + meg +  'Morphed_ico_SemLoc_ica_equalized_Abstract_Mlttap_new_ctf_Coherence_Gamma_150_450_200ms_' + thislabelname
##    print 'gamma'
##    coh_stc.save(data_out)
#	
##############
#
#    print "spectral connectivity ppc"
#    dw_pli=np.zeros((len(labellist),len(labellist),4,2))
#    for actc,tc in enumerate(range(150,251,100)): #time count
#        tminw=tc/1000.
#        tmaxw=(tc+200)/1000.
#        dw_pli[:,:,:,actc], freqs1, times1, n_epochs1, _ = spectral_connectivity(seed_ts, method='ppc', mode='multitaper', sfreq=sfreq, fmin=fmin, fmax=fmax, faverage=True, tmin=tminw, tmax=tmaxw, n_jobs=4)
###spectral_connectivity(data, method='coh', indices=None, sfreq=6.283185307179586, mode='multitaper', fmin=None, fmax=inf, fskip=0, faverage=False, tmin=None, tmax=None, mt_bandwidth=None, mt_adaptive=False, mt_low_bias=True, cwt_frequencies=None, cwt_n_cycles=7, block_size=1000, n_jobs=1, verbose=None)
#    print('ppc fin ')
#    ppc=dw_pli.copy()
#    path_mtx=data_path + meg +'ppc_14ROIs_2_abstract.mat'
#    scio.savemat(path_mtx,{'ppc':ppc})
##
### Generate a SourceEstimate with the Coherence. This is simple since we
### used a single seed. For more than one seeds we would have to split coh.
### Note: We use a hack to save the frequency axis as time
##    print "use a hack to save the frequency axis as time"
##    tmin1 = 0.150#np.mean(freqs[0])
##    tstep1 = 0.100#np.mean(freqs[1]) - tmin
##    coh_stc = mne.SourceEstimate(dw_pli[:,0,:], vertices=stc_coh.vertices, tmin= tmin1, tstep= tstep1, subject='fsaverage')
##    data_out = data_path + meg +  'Morphed_ico_SemLoc_ica_equalized_Abstract_Mlttap_new_ctf_PPC_Theta_150_450_200ms_' + thislabelname
##    print 'theta'
##    coh_stc.save(data_out)
##    
##    coh_stc = mne.SourceEstimate(dw_pli[:,1,:], vertices=stc_coh.vertices, tmin= tmin1, tstep= tstep1, subject='fsaverage')
##    data_out = data_path + meg +  'Morphed_ico_SemLoc_ica_equalized_Abstract_Mlttap_new_ctf_PPC_Alpha_150_450_200ms_' + thislabelname
##    print 'alpha'
##    coh_stc.save(data_out)
##    
##    coh_stc = mne.SourceEstimate(dw_pli[:,2,:], vertices=stc_coh.vertices, tmin= tmin1, tstep= tstep1, subject='fsaverage')
##    data_out = data_path + meg +  'Morphed_ico_SemLoc_ica_equalized_Abstract_Mlttap_new_ctf_PPC_Beta_150_450_200ms_' + thislabelname
##    print 'beta'
##    coh_stc.save(data_out)
##    
##    coh_stc = mne.SourceEstimate(dw_pli[:,3,:], vertices=stc_coh.vertices, tmin= tmin1, tstep= tstep1, subject='fsaverage')
##    data_out = data_path + meg +  'Morphed_ico_SemLoc_ica_equalized_Abstract_Mlttap_new_ctf_PPC_Gamma_150_450_200ms_' + thislabelname
##    print 'gamma'
##    coh_stc.save(data_out)
#
#################
#    print "temp conn mi alpha"
#    
#    MIa=np.zeros((len(labellist),len(labellist),len(stcs),2))
#        
##        for ccc in range(len(stcs)):
##            print ccc		
##            data_in = data_path + meg + 'Morphed_SemLoc_epochs' + str(ccc) + '_icaclean_Abstract_Source_Evoked_m500_700'
##            stc_orig=mne.read_source_estimate(data_in)
##            #stc_orig.resample(100)
##            stcs[ccc]=stc_orig
#        #seed_ts = mne.extract_label_time_course(stcs, label_new, src, mode='max')
#        #comb_ts = zip(seed_ts, stcs)
#        
#        
#    for ccc in range(len(stcs)): #epoch count
#        print ccc		
#        data_in = data_path + meg + 'Morphed_ico_SemLoc_epochs' + str(ccc) + '_icaclean_Abstract_Source_Alpha_m500_700'
#        stc_orig=mne.read_source_estimate(data_in)            
#        seed_ts = mne.extract_label_time_course(stc_orig, labels_new, src, mode='mean_flip')
#            
#        for actc,tc in enumerate(range(150,251,100)): #time count
#            print tc, actc
#            for vc1 in range(MIa.shape[0]):
#                curr_sig=seed_ts[vc1,tc+500:tc+500+200].copy()                
#                for vc2 in range(MIa.shape[1]): #vertex count
#                    if vc1>vc2:
#                        MIa[vc1,vc2,ccc,actc]=mi(curr_sig,seed_ts[vc2,tc+500:tc+500+200],bins=7)#[0,1]                    
#                
#    MIa_final=MIa.mean(axis=2)
#    path_mtx=data_path + meg +'MIa_14ROIs_2_abstract.mat'
#    scio.savemat(path_mtx,{'MIa_final':MIa_final})
#    
#    
#################
#    print "temp conn mi beta"
#    
#    MIb=np.zeros((len(labellist),len(labellist),len(stcs),2))
#        
##        for ccc in range(len(stcs)):
##            print ccc		
##            data_in = data_path + meg + 'Morphed_SemLoc_epochs' + str(ccc) + '_icaclean_Abstract_Source_Evoked_m500_700'
##            stc_orig=mne.read_source_estimate(data_in)
##            #stc_orig.resample(100)
##            stcs[ccc]=stc_orig
#        #seed_ts = mne.extract_label_time_course(stcs, label_new, src, mode='max')
#        #comb_ts = zip(seed_ts, stcs)
#        
#        
#    for ccc in range(len(stcs)): #epoch count
#        print ccc		
#        data_in = data_path + meg + 'Morphed_ico_SemLoc_epochs' + str(ccc) + '_icaclean_Abstract_Source_Beta_m500_700'
#        stc_orig=mne.read_source_estimate(data_in)            
#        seed_ts = mne.extract_label_time_course(stc_orig, labels_new, src, mode='mean_flip')
#            
#        for actc,tc in enumerate(range(150,251,100)): #time count
#            print tc, actc
#            for vc1 in range(MIb.shape[0]):
#                curr_sig=seed_ts[vc1,tc+500:tc+500+200].copy()                
#                for vc2 in range(MIb.shape[1]): #vertex count
#                    if vc1>vc2:
#                        MIb[vc1,vc2,ccc,actc]=mi(curr_sig,seed_ts[vc2,tc+500:tc+500+200],bins=7)#[0,1]                    
#                
#    MIb_final=MIb.mean(axis=2)
#    path_mtx=data_path + meg +'MIb_14ROIs_2_abstract.mat'
#    scio.savemat(path_mtx,{'MIb_final':MIb_final})
#    
#################
#    print "temp conn mi theta"
#    
#    MIt=np.zeros((len(labellist),len(labellist),len(stcs),2))
#        
##        for ccc in range(len(stcs)):
##            print ccc		
##            data_in = data_path + meg + 'Morphed_SemLoc_epochs' + str(ccc) + '_icaclean_Abstract_Source_Evoked_m500_700'
##            stc_orig=mne.read_source_estimate(data_in)
##            #stc_orig.resample(100)
##            stcs[ccc]=stc_orig
#        #seed_ts = mne.extract_label_time_course(stcs, label_new, src, mode='max')
#        #comb_ts = zip(seed_ts, stcs)
#        
#        
#    for ccc in range(len(stcs)): #epoch count
#        print ccc		
#        data_in = data_path + meg + 'Morphed_ico_SemLoc_epochs' + str(ccc) + '_icaclean_Abstract_Source_Theta_m500_700'
#        stc_orig=mne.read_source_estimate(data_in)            
#        seed_ts = mne.extract_label_time_course(stc_orig, labels_new, src, mode='mean_flip')
#            
#        for actc,tc in enumerate(range(150,251,100)): #time count
#            print tc, actc
#            for vc1 in range(MIt.shape[0]):
#                curr_sig=seed_ts[vc1,tc+500:tc+500+200].copy()                
#                for vc2 in range(MIt.shape[1]): #vertex count
#                    if vc1>vc2:
#                        MIt[vc1,vc2,ccc,actc]=mi(curr_sig,seed_ts[vc2,tc+500:tc+500+200],bins=7)#[0,1]                    
#                
#    MIt_final=MIt.mean(axis=2)
#    path_mtx=data_path + meg +'MIt_14ROIs_2_abstract.mat'
#    scio.savemat(path_mtx,{'MIt_final':MIt_final})
#
#################
#    print "temp conn mi gamma"
#    
#    MIg=np.zeros((len(labellist),len(labellist),len(stcs),2))
#        
##        for ccc in range(len(stcs)):
##            print ccc		
##            data_in = data_path + meg + 'Morphed_SemLoc_epochs' + str(ccc) + '_icaclean_Abstract_Source_Evoked_m500_700'
##            stc_orig=mne.read_source_estimate(data_in)
##            #stc_orig.resample(100)
##            stcs[ccc]=stc_orig
#        #seed_ts = mne.extract_label_time_course(stcs, label_new, src, mode='max')
#        #comb_ts = zip(seed_ts, stcs)
#        
#        
#    for ccc in range(len(stcs)): #epoch count
#        print ccc		
#        data_in = data_path + meg + 'Morphed_ico_SemLoc_epochs' + str(ccc) + '_icaclean_Abstract_Source_Gamma_m500_700'
#        stc_orig=mne.read_source_estimate(data_in)            
#        seed_ts = mne.extract_label_time_course(stc_orig, labels_new, src, mode='mean_flip')
#            
#        for actc,tc in enumerate(range(150,251,100)): #time count
#            print tc, actc
#            for vc1 in range(MIg.shape[0]):
#                curr_sig=seed_ts[vc1,tc+500:tc+500+200].copy()                
#                for vc2 in range(MIg.shape[1]): #vertex count
#                    if vc1>vc2:
#                        MIg[vc1,vc2,ccc,actc]=mi(curr_sig,seed_ts[vc2,tc+500:tc+500+200],bins=7)#[0,1]                    
#                
#    MIg_final=MIg.mean(axis=2)
#    path_mtx=data_path + meg +'MIg_14ROIs_2_abstract.mat'
#    scio.savemat(path_mtx,{'MIg_final':MIg_final})