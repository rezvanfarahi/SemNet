# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 15:51:53 2016

@author: rf02
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
#sys.path.append('/imaging/local/software/mne_python/la_v0.9')
sys.path.insert(1,'/imaging/rf02/mne_python_11')
#sys.path.insert(1,'/imaging/local/software/mne_python/la_v0.9full')
#sys.path.insert(1,'/imaging/local/software/mne_python/v0.11')


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
from mne.minimum_norm import (read_inverse_operator, make_inverse_operator, write_inverse_operator, apply_inverse, apply_inverse_epochs)
from mne.minimum_norm.inverse import (prepare_inverse_operator)
import scipy.io
from mne import find_events
from mne.connectivity import seed_target_indices, spectral_connectivity
from mne import read_labels_from_annot
from mne.label import _n_colors as nc
from mne.label import _read_annot as mnera
from mne.label import _write_annot as mnewa
from mne.viz import circular_layout, plot_connectivity_circle
import matplotlib.pyplot as plt
import copy
import pickle
print mne.__path__

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

label_path='/imaging/rf02/TypLexMEG/parcellation/labels_and_CTFs/'
results_path='/imaging/rf02/parcellation/leakage_pattern/'
#labellist_path = [label_path+'rostfront_smallseed-lh.label',label_path+'midtemp_smallseed-lh.label',label_path+'latocci_smallseed-lh.label']#[ 'lh.posterior_fusiform','rh.posterior_fusiform', 'lh.lateraloccipital','rh.lateraloccipital','lh.precentral','rh.precentral', 'lh.ATL_latmed','rh.ATL_latmed']#labellist = [ 'leftATL-lh','rightATL-rh']#, 'medialtempright-rh','medialtempleft-lh']'lh.postcentral','rh.postcentral',
#labellist_path = [label_path+'IFG_smallseed-lh.label',label_path+'IPL_smallseed-lh.label']#[ 'lh.posterior_fusiform','rh.posterior_fusiform', 'lh.lateraloccipital','rh.lateraloccipital','lh.precentral','rh.precentral', 'lh.ATL_latmed','rh.ATL_latmed']#labellist = [ 'leftATL-lh','rightATL-rh']#, 'medialtempright-rh','medialtempleft-lh']'lh.postcentral','rh.postcentral',
labellist_path = [label_path+'LO_smallseed-lh.label',label_path+'RMF_smallseed-lh.label',label_path+'STS_smallseed-lh.label']#[ 'lh.posterior_fusiform','rh.posterior_fusiform', 'lh.lateraloccipital','rh.lateraloccipital','lh.precentral','rh.precentral', 'lh.ATL_latmed','rh.ATL_latmed']#labellist = [ 'leftATL-lh','rightATL-rh']#, 'medialtempright-rh','medialtempleft-lh']'lh.postcentral','rh.postcentral',
labellist=[mne.read_label(label) for label in labellist_path]
for lbl in labellist:
    lbl.values.fill(1.0)
#
#which_parc=0 #names below
#which_conn=1
plt.close("all")

con_names=['COH No Noise No Leakage','COH With Noise No Leakage','COH With Noise With Leakage','imCOH With Noise No Leakage','imCOH With Noise With Leakage']
supertitles=['Desikan-Killiany Atlas (68 ROIs)','Modified Desikan-Killiany Atlas (74 ROIs)','Destrieux Atlas (148 ROIs)','Modified Destrieux Atlas (74 ROIs)','Region Growing (70 ROIs)']

for which_parc in np.array([2]):#range(5):#np.array([0]):# in range(1):
    print which_parc
    fig = plt.figure(num=which_parc, figsize=(20, 6), facecolor='black')
    fig.suptitle(supertitles[which_parc],color='white',fontsize=18)
    for which_conn in range(5):#np.array([0,2]):#
        print which_conn
        parc_name=['aparc',
        'fsaverage_finalpostparc_smooth_stddiff_split_merge_halfmax_86_badremoved_74',
        'aparc.a2009s',
        'fsaverage_finalpostparc2009_smooth_stddiff_split_merge_halfmax_98_badremoved_74_cmatch',
        'fsaverage_RG_finalpostparc_smooth_stddiff_split_merge_halfmax_82_badremoved_70_cmatch'
        ]##'fsaverage_finalpostparc2009_smooth_stddiff_split_merge_halfmax_98_badremoved_74_cmatch'#'fsaverage_RG_finalpostparc_smooth_stddiff_split_merge_halfmax_82_badremoved_70_cmatch'
        #
        labels_pre = mne.read_labels_from_annot('fsaverage', parc=parc_name[which_parc], subjects_dir=subjects_dir)
        #labels_pre.remove(labels_pre[-1])
        #    labels_pre.remove(labels_pre[-1])
        Rmat_names=['ctf_ypos_aparc.mat',
        'ctf_ypos_finalpostparc_smooth_stddiff_split_merge_halfmax_86_badremoved_74.mat',
        'ctf_ypos_anat2009.mat',
        'ctf_ypos_finalpostparc2009_smooth_stddiff_split_merge_halfmax_98_badremoved_74_cmatch.mat',
        'ctf_ypos_RG_finalpostparc_smooth_stddiff_split_merge_halfmax_82_badremoved_70_cmatch.mat'
        ]
        path_Rmat=results_path+Rmat_names[which_parc]
        ltc=scipy.io.loadmat(path_Rmat,mat_dtype=True); ctf_to_ypos=np.squeeze(ltc['ctf_to_ypos'])
        
        ypos_names=['ypos_aparc.mat',
        'ypos_aparcmod_badremoved_74.mat',
        'ypos_aparc2009.mat',
        'ypos_aparc2009mod_badremoved_74.mat',
        'ypos_RG_badremoved_70.mat'
        ]
        path_ypos=results_path+ypos_names[which_parc]
        ltc=scipy.io.loadmat(path_ypos,mat_dtype=True); ypos_sort=np.squeeze(ltc['ypos_sort'])
        labels=[labels_pre[jj] for jj in ypos_sort]
        inname=['aparc','aparc_mod','aparc2009','aparc2009_mod','RG']
        fname_a=['nonoise','_ideal_coh','_coh','_ideal_imcoh','_imcoh']#%'aparc_coh_mod.mat';
        if which_conn==0:
            fname=label_path+inname[which_parc]+fname_a[which_conn+1]+'prob_mat_hubs.mat';
            ltc=scipy.io.loadmat(fname,mat_dtype=True); prob_mat=np.squeeze(ltc['probmat_hubs'])
            prob_mat[prob_mat<0.5]=0
        else:
            fname=label_path+inname[which_parc]+fname_a[which_conn]+'prob_mat_hubs.mat';
            ltc=scipy.io.loadmat(fname,mat_dtype=True); prob_mat=np.squeeze(ltc['probmat_hubs'])
        
        label_names = [label.name for label in labels]
        label_colors1 = [label.color for label in labels]
        label_colors=[tuple((label_colors1[iii][0],label_colors1[iii][1],label_colors1[iii][2],1.0)) for iii in range(len(label_colors1))]
        lh_labels = [name for name in label_names if name.endswith('lh')]
        
        # Get the y-location of the label
        label_ypos = list()
        for name in lh_labels:
            idx = label_names.index(name)
            ypos = np.mean(labels[idx].pos[:, 1])
            label_ypos.append(ypos)
        
        # Reorder the labels based on their location
        lh_labels = [label for (yp, label) in sorted(zip(label_ypos, lh_labels))]
        
        # For the right hemi
        rh_labels = [label[:-2] + 'rh' for label in lh_labels]
        
        # Save the plot order and create a circular layout
        node_order = list()
        node_order.extend(lh_labels[::-1])  # reverse the order
        node_order.extend(rh_labels)
        
        node_angles = circular_layout(label_names, node_order, start_pos=90,
                                      group_boundaries=[0, len(label_names) / 2])
        
        # Plot the graph using node colors from the FreeSurfer parcellation. We only
        # show the 300 strongest connections.
        plot_connectivity_circle(np.abs(prob_mat), range(len(label_names)),colormap='hot',n_lines=300,
                                 node_angles=node_angles, node_colors=label_colors,vmin=0.,vmax=1.,colorbar=True, fontsize_title=14.,fontsize_names=6.,fontsize_colorbar=12.,colorbar_size=0.4,
                                 fig=fig, subplot=(1, 5, which_conn+1),padding=1.0, title=con_names[which_conn])
        
        
        #        outname=['aparc','aparc_mod','aparc2009','aparc2009_mod','RG']
        
        plt.show(block=False)
    outname=label_path+inname[which_parc]+'allconn_new2_cbar.jpg'                
    plt.savefig(outname,bbox_inches='tight',facecolor=fig.get_facecolor(),transparent=True)

