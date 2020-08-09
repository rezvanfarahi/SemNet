"""
Created on Tue Jul 11 06:26:19 2017

@author: rf02
"""

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
sys.path.insert(1,'/imaging/rf02/TypLexMEG/mne_python0.9')

# for qsub
# (add ! here if needed) 
sys.path.append('/imaging/local/software/anaconda/latest/x86_64/bin/python')
sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
sys.path.insert(1,'/imaging/local/software/anaconda/2.4.1/3/lib/python3.5/site-packages/scipy')
###

sys.path.insert(1,'/imaging/local/software/python_packages/scikit-learn/v0.13.1/lib.linux-x86_64-2.7/sklearn/externals')
#import joblib
###


import os
import numpy as np
import mne
import scipy
from scipy.stats.mstats import zscore
import scipy.io as scio
#import sklearn
from mne.io import Raw
from mne.epochs import equalize_epoch_counts
from mne.minimum_norm import (read_inverse_operator, psf_ctf_28Oct15, apply_inverse)
from mne.stats import f_threshold_mway_rm, f_mway_rm
from mne import (io, spatial_tris_connectivity, grade_to_tris)
from mne.epochs import equalize_epoch_counts
from mne.stats import (permutation_cluster_test,summarize_clusters_stc)
from mne.minimum_norm import apply_inverse, read_inverse_operator
from mne.stats import f_threshold_mway_rm, f_mway_rm, fdr_correction, spatio_temporal_cluster_test
from mne.datasets import sample
from mne.viz import mne_analyze_colormap
import copy

###############################################################################
#out_path = '/home/rf02/rezvan/test1/step_by_step/dcm/' # root directory for your MEG data
out_path = '/imaging/rf02/Semnet/semnet4semloc/'

subjects_dir = '/imaging/rf02/Semnet/'    # where your MRI subdirectories are
# where event-files aremean
data_path = '/imaging/rf02/Semnet/'    # where event files are

label_path = '/home/rf02/rezvan/TypLexMEG/createdlabels_SL/'

inv_path = '/imaging/rf02/Semnet/'

# get indices for subjects to be processed from command line input
# 
print (sys.argv)
subject_inds = []
for ss in sys.argv[1:]:
   subject_inds.append( int( ss ) )

subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]
print ("subject_inds:")
print (subject_inds)
print ("No rejection")

list_all =  ['/meg16_0030/160216/', #0
            '/meg16_0032/160218/', #1
            '/meg16_0034/160219/', #3
            '/meg16_0035/160222/', #4
            '/meg16_0042/160229/', #7
            '/meg16_0045/160303/', #8
            '/meg16_0052/160310/', #10
            '/meg16_0056/160314/',#11
            '/meg16_0069/160405/',#12 
            '/meg16_0070/160407/', #13
            '/meg16_0072/160408/', #15
            '/meg16_0073/160411/', #16
            '/meg16_0075/160411/', #17
            '/meg16_0078/160414/', #18
            '/meg16_0082/160418/', #19
            '/meg16_0086/160422/', #20
            '/meg16_0097/160512/', #21 
            '/meg16_0122/160707/', #22 
            '/meg16_0125/160712/', #24 
            ]



# subjects names used for MRI data
subjects=  ['MRI_meg16_0030' ,#0
            'MRI_meg16_0032' ,#1
            'MRI_meg16_0034' ,#2
            'MRI_meg16_0035' ,#3
            'MRI_meg16_0042' ,#4
            'MRI_meg16_0045' ,#5
            'MRI_meg16_0052' ,#6
            'MRI_meg16_0056' ,#7
    	       'MRI_meg16_0069' ,#8
 	       'MRI_meg16_0070' ,#9
            'MRI_meg16_0072' ,#10
            'MRI_meg16_0073' ,#11
            'MRI_meg16_0075' ,#12
            'MRI_meg16_0078' ,#13
            'MRI_meg16_0082' ,#14
            'MRI_meg16_0086' ,#15
            'MRI_meg16_0097' ,#16
            'MRI_meg16_0122' ,#17
            'MRI_meg16_0125' ,#18
            ]

ll = []
for ss in subject_inds:
 ll.append(list_all[ss])

print ("ll:")
print (ll)
labelss = mne.read_labels_from_annot(subject='fsaverage', parc='rez_aparc_25oct16', subjects_dir=data_path)

label_names=[label.name for label in labelss]
label_index=np.array([label_names.index('ATL_latmed-lh'),label_names.index('posterior_inferiorfrontal-lh'),label_names.index('middle_temporal-lh'),label_names.index('AG_new-lh'),label_names.index('vWFA2-lh')])#label_names.index('ATL_latmed-lh'),label_names.index('supramarginal-lh'),
lfname=data_path+'fsaverage/label/ATL_manual_prefinal2-lh.label'
labelatllh=mne.read_label(lfname,subject='fsaverage')

labellist=[labelss[ii] for ii in label_index]
#labellist[0]=labelatllh.copy()
labellist=[labelatllh,labelss[label_index[1]],labelss[label_index[2]],labelss[label_index[3]]]
for lbl in labellist:
    lbl.values.fill(1.0)


srcin='/imaging/rf02/TypLexMEG/fsaverage/bem/fsaverage_dist-ico-5-src.fif'
src_avg = mne.read_source_spaces(srcin)
fname = data_path + list_all[0] + 'firstMorphed_ico_oldreg_LD_SL_1_48ica_Concrete_Source_Signed_Evoked_m300_600'
this_stc = mne.read_source_estimate(fname)  
tmin1=-300
tstep1=1
refstc_data=np.ones((this_stc.data.shape[0],1))
vertices_to = [np.arange(10242), np.arange(10242)]
refstc = mne.SourceEstimate(refstc_data, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
#thisstc_label=mne.stc_to_label(thisstc_tolabel, src_avg)

ntimes=refstc.data.shape[1]
for this_labelc, this_label in enumerate(labellist): 
    label_origdata=refstc.in_label(this_label).data
    label_verts=refstc.in_label(this_label).lh_vertno
    refstc_data[label_verts,0]=0
    thisstc_tolabel = mne.SourceEstimate(refstc_data, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
    thisstc_label=mne.stc_to_label(thisstc_tolabel, src_avg,smooth=True)[0]
    label_path=out_path+'/'+'mask_labels_ATL_IFG_MTG_AG-lh'
    thisstc_label.save(label_path)

#    vertices_test = [np.arange(10242), np.arange(10242)]
#            np.intersect1d(label_verts,vertices_to[0])
#    vert50=np.sort(Matx[label_verts,this_labelc])[-1]
#    winner_verts=label_verts[np.where(Matx[label_verts,this_labelc]>=vert50)[0]]#label_verts[np.where(Matx[label_verts,this_labelc]>=np.max(Matx[label_verts,this_labelc])/2.)[0]]#np.where(Matx[:,this_labelc]>=np.max(Matx[:,this_labelc])/2.)[0]#
#    #            print winner_verts
#    thisscale=np.hstack((mne.label.label_sign_flip(thisstc_label[0],src_avg),mne.label.label_sign_flip(thisstc_label[1],src_avg)))[winner_verts]#mne.label.label_sign_flip(this_label,src_avg)[winner_verts]#1#Matx[winner_verts,this_labelc]/np.max(Matx[winner_verts,this_labelc])
#    label_predata=this_stc.data[winner_verts,:]
#    #            label_data1=np.tile(thisscale,[ntimes,1]).transpose()*label_predata
#    label_data1=thisscale[:,np.newaxis]*label_predata
#    active_ones=np.where(np.mean(np.abs(label_data1[:,sonset:]),1)>=np.min(np.min(np.abs(label_data1[:,sonset:]),1)))[0]
#    print (len(active_ones))
#    label_data=label_data1[active_ones,:]
#    ini_tc=np.mean(label_data[:,:], axis=0)
#    thisstc=np.zeros((this_stc.data.shape[0],1))
#    thisstc[label_verts,0]=1