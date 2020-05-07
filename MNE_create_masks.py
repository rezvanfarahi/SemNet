# -*- coding: utf-8 -*-
"""
Created on Tue Jan  3 19:15:24 2017

@author: rf02
"""

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
#sys.path.insert(1,'/imaging/local/software/mne_python/latest_v0.9')
sys.path.append('/imaging/local/software/python_packages/scikit-learn/v0.16.1')
sys.path.insert(1,'/imaging/local/software/python_packages/scikit-learn/v0.13.1/lib.linux-x86_64-2.7/sklearn/externals')
import joblib
###


import os
import numpy as np
import mne
import scipy
from scipy.stats.mstats import zscore
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
import matplotlib.pyplot as plt
main_path = '/imaging/rf02/Semnet/'
data_path = '/imaging/rf02/Semnet'	# where subdirs for MEG data are
event_path = '/imaging/rf02/Semnet'
os.chdir(data_path)
label_path = '/imaging/rf02/Semnet/masks/'#'/imaging/rf02/TypLexMEG/fsaverage/label/'
inv_path = '/imaging/rf02/Semnet/'
subjects_dir = data_path
print sys.argv
subject_inds = []
for ss in sys.argv[1:]:
   subject_inds.append( int( ss ) )

subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]
print "subject_inds:"
print subject_inds
print "No rejection"

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

print "ll:"
print ll

labellist_path = [label_path+'pfusiform_L-lh.label',label_path+'pfusiform_R-rh.label',label_path+'hand_final_L-lh.label',label_path+'hand_final_R-rh.label']#[ 'lh.posterior_fusiform','rh.posterior_fusiform', 'lh.lateraloccipital','rh.lateraloccipital','lh.precentral','rh.precentral', 'lh.ATL_latmed','rh.ATL_latmed']#labellist = [ 'leftATL-lh','rightATL-rh']#, 'medialtempright-rh','medialtempleft-lh']'lh.postcentral','rh.postcentral',
labellist1=[mne.read_label(label) for label in labellist_path]
for thisl in labellist1:
    thisl.values.fill(1.0)
    thisl.subject='fsaverage'
fsaverage_path='/imaging/rf02/TypLexMEG'
labelss1 = mne.read_labels_from_annot(subject='fsaverage', parc='aparc', subjects_dir=fsaverage_path)
#labelss.remove(labelss[-3])
#labelss.remove(labelss[-2])
labelss1.remove(labelss1[-1])
for lbl in labelss1:
    lbl.values.fill(1.0)
labelss=labelss1+labellist1
label_names=[label.name for label in labelss]
MT=np.ones((20484,1))
vertices_to = [np.arange(10242), np.arange(10242)]
this_stc = mne.SourceEstimate(MT, vertices=vertices_to,tmin=1e-3 * 0, tstep=1e-3 * 1, subject='fsaverage')
srcin='/imaging/rf02/TypLexMEG/fsaverage/bem/fsaverage_dist-ico-5-src.fif'
#src_avg2=src_avg.copy()
src_avg = mne.read_source_spaces(srcin)

#label_index_vishand=np.array([label_names.index('precentral-lh'),label_names.index('postcentral-lh'),label_names.index('parsopercularis-lh'),label_names.index('parsorbitalis-lh'),label_names.index('parstriangularis-lh'),label_names.index('lateraloccipital-lh'),label_names.index('lingual-lh'),label_names.index('fusiform-lh'),label_names.index('pericalcarine-lh'),label_names.index('cuneus-lh'),label_names.index('inferiortemporal-lh')])
#mask_vishand_L=labelss[label_index_vishand[0]]
#for llcnt in label_index_vishand[1:]:
#    mask_vishand_L=mask_vishand_L+labelss[llcnt]
#thismask_verts=this_stc.in_label(mask_vishand_L)
#MT=np.ones((20484,1))
#MT[thismask_verts.lh_vertno]=0
#thismask_stc = mne.SourceEstimate(MT, vertices=vertices_to,tmin=1e-3 * 0, tstep=1e-3 * 1, subject='fsaverage')
#excludemask=mne.stc_to_label(thismask_stc,src=src_avg, smooth=True,connected=False, subjects_dir=fsaverage_path)
#excludemaskl=excludemask[0]
#excludemaskl.save('/imaging/rf02/Semnet/masks/mask_vishand_L-lh.label')
#    
#label_index_vishear=np.array([label_names.index('transversetemporal-lh'),label_names.index('superiortemporal-lh'),label_names.index('bankssts-lh'),label_names.index('insula-lh'),label_names.index('lateraloccipital-lh'),label_names.index('lingual-lh'),label_names.index('fusiform-lh'),label_names.index('pericalcarine-lh'),label_names.index('cuneus-lh'),label_names.index('inferiortemporal-lh')])
#mask_vishear_L=labelss[label_index_vishear[0]]
#for llcnt in label_index_vishear[1:]:
#    mask_vishear_L=mask_vishear_L+labelss[llcnt]
#thismask_verts=this_stc.in_label(mask_vishear_L)
#MT=np.ones((20484,1))
#MT[thismask_verts.lh_vertno]=0
#thismask_stc = mne.SourceEstimate(MT, vertices=vertices_to,tmin=1e-3 * 0, tstep=1e-3 * 1, subject='fsaverage')
#excludemask=mne.stc_to_label(thismask_stc,src=src_avg, smooth=True,connected=False, subjects_dir=fsaverage_path)
#excludemaskl=excludemask[0]
#excludemaskl.save('/imaging/rf02/Semnet/masks/mask_vishear_L-lh.label')
#    
#label_index_handhear=np.array([label_names.index('transversetemporal-lh'),label_names.index('superiortemporal-lh'),label_names.index('bankssts-lh'),label_names.index('insula-lh'),label_names.index('precentral-lh'),label_names.index('postcentral-lh'),label_names.index('parsopercularis-lh'),label_names.index('parsorbitalis-lh'),label_names.index('parstriangularis-lh')])
#mask_handhear_L=labelss[label_index_handhear[0]]
#for llcnt in label_index_handhear[1:]:
#    mask_handhear_L=mask_handhear_L+labelss[llcnt]
#thismask_verts=this_stc.in_label(mask_handhear_L)
#MT=np.ones((20484,1))
#MT[thismask_verts.lh_vertno]=0
#thismask_stc = mne.SourceEstimate(MT, vertices=vertices_to,tmin=1e-3 * 0, tstep=1e-3 * 1, subject='fsaverage')
#excludemask=mne.stc_to_label(thismask_stc,src=src_avg, smooth=True,connected=False, subjects_dir=fsaverage_path)
#excludemaskl=excludemask[0]
#excludemaskl.save('/imaging/rf02/Semnet/masks/mask_handhear_L-lh.label')    
#
#label_index_vishand=np.array([label_names.index('precentral-rh'),label_names.index('postcentral-rh'),label_names.index('parsopercularis-rh'),label_names.index('parsorbitalis-rh'),label_names.index('parstriangularis-rh'),label_names.index('lateraloccipital-rh'),label_names.index('lingual-rh'),label_names.index('fusiform-rh'),label_names.index('pericalcarine-rh'),label_names.index('cuneus-rh'),label_names.index('inferiortemporal-rh')])
#mask_vishand_R=labelss[label_index_vishand[0]]
#for llcnt in label_index_vishand[1:]:
#    mask_vishand_R=mask_vishand_R+labelss[llcnt]
#thismask_verts=this_stc.in_label(mask_vishand_R)
#MT=np.ones((20484,1))
#MT[thismask_verts.rh_vertno+10242]=0
#thismask_stc = mne.SourceEstimate(MT, vertices=vertices_to,tmin=1e-3 * 0, tstep=1e-3 * 1, subject='fsaverage')
#excludemask=mne.stc_to_label(thismask_stc,src=src_avg, smooth=True,connected=False, subjects_dir=fsaverage_path)
#excludemaskr=excludemask[1]
#excludemaskr.save('/imaging/rf02/Semnet/masks/mask_vishand_R-rh.label')
#
#    
#label_index_vishear=np.array([label_names.index('transversetemporal-rh'),label_names.index('superiortemporal-rh'),label_names.index('bankssts-rh'),label_names.index('insula-rh'),label_names.index('lateraloccipital-rh'),label_names.index('lingual-rh'),label_names.index('fusiform-rh'),label_names.index('pericalcarine-rh'),label_names.index('cuneus-rh'),label_names.index('inferiortemporal-rh')])
#mask_vishear_R=labelss[label_index_vishear[0]]
#for llcnt in label_index_vishear[1:]:
#    mask_vishear_R=mask_vishear_R+labelss[llcnt]
#thismask_verts=this_stc.in_label(mask_vishear_R)
#MT=np.ones((20484,1))
#MT[thismask_verts.rh_vertno+10242]=0
#thismask_stc = mne.SourceEstimate(MT, vertices=vertices_to,tmin=1e-3 * 0, tstep=1e-3 * 1, subject='fsaverage')
#excludemask=mne.stc_to_label(thismask_stc,src=src_avg, smooth=True,connected=False, subjects_dir=fsaverage_path)
#excludemaskr=excludemask[1]
#excludemaskr.save('/imaging/rf02/Semnet/masks/mask_vishear_R-rh.label')
#
#    
#label_index_handhear=np.array([label_names.index('transversetemporal-rh'),label_names.index('superiortemporal-rh'),label_names.index('bankssts-rh'),label_names.index('insula-rh'),label_names.index('precentral-rh'),label_names.index('postcentral-rh'),label_names.index('parsopercularis-rh'),label_names.index('parsorbitalis-rh'),label_names.index('parstriangularis-rh')])
#mask_handhear_R=labelss[label_index_handhear[0]]
#for llcnt in label_index_handhear[1:]:
#    mask_handhear_R=mask_handhear_R+labelss[llcnt]
#thismask_verts=this_stc.in_label(mask_handhear_R)
#MT=np.ones((20484,1))
#MT[thismask_verts.rh_vertno+10242]=0
#thismask_stc = mne.SourceEstimate(MT, vertices=vertices_to,tmin=1e-3 * 0, tstep=1e-3 * 1, subject='fsaverage')
#excludemask=mne.stc_to_label(thismask_stc,src=src_avg, smooth=True,connected=False, subjects_dir=fsaverage_path)
#excludemaskr=excludemask[1]
#excludemaskr.save('/imaging/rf02/Semnet/masks/mask_handhear_R-rh.label')    




#####################################UNCOMMENT FROM HERE
#label_index_vishand=np.array([label_names.index('hand_final_L-lh'),label_names.index('parsopercularis-lh'),label_names.index('parstriangularis-lh'),label_names.index('lateraloccipital-lh'),label_names.index('pfusiform_L-lh')])
#MT=np.ones((20484,1))
#for llcnt in label_index_vishand:
#    thismask_verts=this_stc.in_label(labelss[llcnt])
#    MT[thismask_verts.lh_vertno]=0
#thismask_stc = mne.SourceEstimate(MT, vertices=vertices_to,tmin=1e-3 * 0, tstep=1e-3 * 1, subject='fsaverage')
#excludemask=mne.stc_to_label(thismask_stc,src=src_avg, smooth=True,connected=False, subjects_dir=fsaverage_path)
#excludemaskl=excludemask[0]
#excludemaskl.save('/imaging/rf02/Semnet/masks/mask_vishand_L-lh.label')

label_index_vishand=np.array([label_names.index('lateralorbitofrontal-lh'),label_names.index('parsorbitalis-lh'),label_names.index('parsopercularis-lh'),label_names.index('parstriangularis-lh'),label_names.index('precentral-lh'),label_names.index('postcentral-lh'),label_names.index('paracentral-lh'),label_names.index('supramarginal-lh'),label_names.index('inferiorparietal-lh'),label_names.index('lateraloccipital-lh'),label_names.index('fusiform-lh'),label_names.index('inferiortemporal-lh'),label_names.index('middletemporal-lh'),label_names.index('superiortemporal-lh'),label_names.index('transversetemporal-lh'),label_names.index('bankssts-lh'),label_names.index('temporalpole-lh'),label_names.index('entorhinal-lh'),label_names.index('parahippocampal-lh'),label_names.index('lingual-lh'),label_names.index('pericalcarine-lh'),label_names.index('cuneus-lh'),label_names.index('medialorbitofrontal-lh')])
MT=np.ones((20484,1))
for llcnt in label_index_vishand:
    thismask_verts=this_stc.in_label(labelss[llcnt])
    MT[thismask_verts.lh_vertno]=0
thismask_stc = mne.SourceEstimate(MT, vertices=vertices_to,tmin=1e-3 * 0, tstep=1e-3 * 1, subject='fsaverage')
excludemask=mne.stc_to_label(thismask_stc,src=src_avg, smooth=True,connected=False, subjects_dir=fsaverage_path)
excludemaskl=excludemask[0]
excludemaskl.save('/imaging/rf02/Semnet/masks/mask_vishand_lang_L-lh.label')
excludemaskr=excludemask[1] 
excludemaskr.save('/imaging/rf02/Semnet/masks/mask_vishand_lang_R-rh.label')
#label_index_vishear=np.array([label_names.index('transversetemporal-lh'),label_names.index('superiortemporal-lh'),label_names.index('bankssts-lh'),label_names.index('insula-lh'),label_names.index('lateraloccipital-lh'),label_names.index('pfusiform_L-lh')])
#mask_vishear_L=labelss[label_index_vishear[0]]
#MT=np.ones((20484,1))
#for llcnt in label_index_vishear:
#    thismask_verts=this_stc.in_label(labelss[llcnt])
#    MT[thismask_verts.lh_vertno]=0
#thismask_stc = mne.SourceEstimate(MT, vertices=vertices_to,tmin=1e-3 * 0, tstep=1e-3 * 1, subject='fsaverage')
#excludemask=mne.stc_to_label(thismask_stc,src=src_avg, smooth=True,connected=False, subjects_dir=fsaverage_path)
#excludemaskl=excludemask[0]
#excludemaskl.save('/imaging/rf02/Semnet/masks/mask_vishear_L-lh.label')
#    
#label_index_handhear=np.array([label_names.index('transversetemporal-lh'),label_names.index('superiortemporal-lh'),label_names.index('bankssts-lh'),label_names.index('insula-lh'),label_names.index('hand_final_L-lh'),label_names.index('parsopercularis-lh'),label_names.index('parstriangularis-lh')])
#MT=np.ones((20484,1))
#for llcnt in label_index_handhear:
#    thismask_verts=this_stc.in_label(labelss[llcnt])
#    MT[thismask_verts.lh_vertno]=0
#thismask_stc = mne.SourceEstimate(MT, vertices=vertices_to,tmin=1e-3 * 0, tstep=1e-3 * 1, subject='fsaverage')
#excludemask=mne.stc_to_label(thismask_stc,src=src_avg, smooth=True,connected=False, subjects_dir=fsaverage_path)
#excludemaskl=excludemask[0]
#excludemaskl.save('/imaging/rf02/Semnet/masks/mask_handhear_L-lh.label')    
#
#label_index_vishand=np.array([label_names.index('hand_final_R-rh'),label_names.index('parsopercularis-rh'),label_names.index('parstriangularis-rh'),label_names.index('lateraloccipital-rh'),label_names.index('pfusiform_R-rh')])
#MT=np.ones((20484,1))
#for llcnt in label_index_vishand:
#    thismask_verts=this_stc.in_label(labelss[llcnt])
#    MT[thismask_verts.rh_vertno+10242]=0
#thismask_stc = mne.SourceEstimate(MT, vertices=vertices_to,tmin=1e-3 * 0, tstep=1e-3 * 1, subject='fsaverage')
#excludemask=mne.stc_to_label(thismask_stc,src=src_avg, smooth=True,connected=False, subjects_dir=fsaverage_path)
#excludemaskr=excludemask[1]
#excludemaskr.save('/imaging/rf02/Semnet/masks/mask_vishand_R-rh.label')
#
#    
#label_index_vishear=np.array([label_names.index('transversetemporal-rh'),label_names.index('superiortemporal-rh'),label_names.index('bankssts-rh'),label_names.index('insula-rh'),label_names.index('lateraloccipital-rh'),label_names.index('pfusiform_R-rh')])
#MT=np.ones((20484,1))
#for llcnt in label_index_vishear:
#    thismask_verts=this_stc.in_label(labelss[llcnt])
#    MT[thismask_verts.rh_vertno+10242]=0
#thismask_stc = mne.SourceEstimate(MT, vertices=vertices_to,tmin=1e-3 * 0, tstep=1e-3 * 1, subject='fsaverage')
#excludemask=mne.stc_to_label(thismask_stc,src=src_avg, smooth=True,connected=False, subjects_dir=fsaverage_path)
#excludemaskr=excludemask[1]
#excludemaskr.save('/imaging/rf02/Semnet/masks/mask_vishear_R-rh.label')
#
#    
#label_index_handhear=np.array([label_names.index('transversetemporal-rh'),label_names.index('superiortemporal-rh'),label_names.index('bankssts-rh'),label_names.index('insula-rh'),label_names.index('hand_final_R-rh'),label_names.index('parsopercularis-rh'),label_names.index('parstriangularis-rh')])
#MT=np.ones((20484,1))
#for llcnt in label_index_handhear:
#    thismask_verts=this_stc.in_label(labelss[llcnt])
#    MT[thismask_verts.rh_vertno+10242]=0
#thismask_stc = mne.SourceEstimate(MT, vertices=vertices_to,tmin=1e-3 * 0, tstep=1e-3 * 1, subject='fsaverage')
#excludemask=mne.stc_to_label(thismask_stc,src=src_avg, smooth=True,connected=False, subjects_dir=fsaverage_path)
#excludemaskr=excludemask[1]
#excludemaskr.save('/imaging/rf02/Semnet/masks/mask_handhear_R-rh.label')    