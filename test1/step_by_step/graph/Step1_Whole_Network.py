"""
=========================================================================
Compute source space connectivity and visualize it using a circular graph
=========================================================================

This example computes the all-to-all connectivity between 68 regions in
source space based on dSPM inverse solutions and a FreeSurfer cortical
parcellation. The connectivity is visualized using a circular graph which
is ordered based on the locations of the regions.
"""

# Authors: Martin Luessi <mluessi@nmr.mgh.harvard.edu>
#          Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
#          Nicolas P. Rougier (graph code borrowed from his matplotlib gallery)
#
# License: BSD (3-clause)

print(__doc__)
import sys

#sys.path.append('/imaging/local/software/python_packages/scikit-learn/v0.14.1')
sys.path.append('/imaging/local/software/mne_python/latest_v0.9')

# for qsub
# (add ! here if needed) 
sys.path.append('/imaging/local/software/EPD/latest/x86_64/bin/python')
#sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
###

import numpy as np
import mne
from mne.io import Raw
from mne.minimum_norm import (read_inverse_operator, apply_inverse, apply_inverse_epochs)
import scipy.io
from mne import find_events
from mne.connectivity import (spectral_connectivity, phase_slope_index)
import numpy as np
from mne import read_labels_from_annot
from mne.viz import circular_layout, plot_connectivity_circle



data_path = '/imaging/rf02/TypLexMEG/' # root directory for your MEG data
subjects_dir = '/imaging/rf02/TypLexMEG/'    # where your MRI subdirectories are
out_path='/imaging/rf02/TypLexMEG/icaanalysis_results/jpg/time_resolved/ROI_connectivity/graph/'

event_path = event_path = '/group/erp/data/olaf.hauk/MEG/TypLexM/data/MF22_new/'   # where event files are

label_path = '/imaging/rf02/TypLexMEG/fsaverage/label/'

inv_path = '/imaging/rf02/TypLexMEG/'
inv_fname = 'InverseOperator_fft_1_48_clean_ica_EMEG-inv.fif'

# get indices for subjects to be processed from command line input
# 
print sys.argv
subject_inds = []
for ss in sys.argv[1:]:
   subject_inds.append( int( ss ) )



subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]#0,1,2,3,4,5,6,7,8]#,
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
Dmat_cncrt=np.zeros((46,46,4,752,17))
Dmat_abs=np.zeros((46,46,4,752,17))
Dmat_sub=np.zeros((46,46,4,752,17))
Dmat_both=np.zeros((46,46,4,752,17))
Dmat_b=np.zeros((46,46,4,12))
Dmat_br=np.zeros((46,46,4,12,17))#int(Dmat_both.shape[3]/5),17))


for cnt1,meg in enumerate(ll):
    print cnt1
    
    path_c=data_path+ meg+'thresholded_15_pca_ppc_concrete.mat'
    ltc=scipy.io.loadmat(path_c,mat_dtype=True); thresh_c=ltc['thresh_c']; 
    #scipy.io.savemat(path_c,{'label_ts_concrete':label_ts_concrete})
    path_a=data_path+ meg+'thresholded_15_pca_ppc_abstract.mat'
    #scipy.io.savemat(path_a,{'label_ts_abstract':label_ts_abstract})
    lta=scipy.io.loadmat(path_a,mat_dtype=True); thresh_a=lta['thresh_a']; 
    path_b=data_path+ meg+'thresholded_15_pca_ppc_both.mat'
    ltb=scipy.io.loadmat(path_b,mat_dtype=True); thresh_b=ltb['thresh_b']; 
    #scipy.io.savemat(path_c,{'label_ts_concrete':label_ts_concrete})
    print "lbl ts"
    Dmat_cncrt[:,:,:,:,cnt1]=thresh_c
    Dmat_abs[:,:,:,:,cnt1]=thresh_a
    Dmat_sub[:,:,:,:,cnt1]=thresh_c-thresh_a
    Dmat_both[:,:,:,:,cnt1]=thresh_b
labels = mne.read_labels_from_annot(subject='fsaverage', parc='rezvan_parc_mne_ctf_mn', subjects_dir=subjects_dir)
labels=labels[:-2]
label_colors = [label.color for label in labels]

path_cmap=data_path+ 'icaanalysis_results/rez_cmap.mat'
cmapp=scipy.io.loadmat(path_cmap,mat_dtype=True); cm=cmapp['rez_cmap']
import matplotlib as mpl
cmap=mpl.colors.ListedColormap(cm/255)
# Now, we visualize the connectivity using a circular graph layout

# First, we reorder the labels based on their location in the left hemi
labels[0].name=labels[0].name.replace('ATL','anterior_temporal')
labels[1].name=labels[1].name.replace('ATL','anterior_temporal')
#labels[1].name=labels[1].name.replace('anterior_temporal','anterior_temporal'+str(5))
label_names = [label.name for label in labels]

lh_labels = [name for name in label_names if name.endswith('lh')]

# Get the y-location of the label
label_ypos = list()
for name in lh_labels:
    idx = label_names.index(name)
    ypos = np.mean(labels[idx].pos[:, 1])
    label_ypos.append(ypos)

# Reorder the labels based on their location
lh_labels = [label for (ypos, label) in sorted(zip(label_ypos, lh_labels))]

# For the right hemi
rh_labels = [label[:-2] + 'rh' for label in lh_labels]

# Save the plot order and create a circular layout
node_order = list()
node_order.extend(lh_labels[::-1])  # reverse the order
node_order.extend(rh_labels)

node_angles = circular_layout(label_names, node_order, start_pos=90, group_boundaries=[0, len(label_names) / 2])
#b2=5
#for cc in range(int(Dmat_both.shape[3]/5)):
#    print cc
#    b1=b2-5; b2=b1+10
#    for cntt1 in range(Dmat_both.shape[0]):
#        for cntt2 in range(Dmat_both.shape[1]):
#            for cntt3 in range(Dmat_both.shape[2]):
#                for cntt4 in range(Dmat_both.shape[4]):
#                    if np.sum(Dmat_both[cntt1,cntt2,cntt3,b1:b2,cntt4]>0)>5:
#                        Dmat_br[cntt1,cntt2,cntt3,cc,cntt4]=1
#                    else:
#                        Dmat_br[cntt1,cntt2,cntt3,cc,cntt4]=0

b2=150
for cc in range(12):
    print cc
    b1=b2-50; b2=b1+100
    for cntt1 in range(Dmat_both.shape[0]):
        for cntt2 in range(Dmat_both.shape[1]):
            for cntt3 in range(Dmat_both.shape[2]):
                for cntt4 in range(Dmat_both.shape[4]):
                    if np.sum(Dmat_both[cntt1,cntt2,cntt3,b1:b2,cntt4]>0)>50:
                        Dmat_br[cntt1,cntt2,cntt3,cc,cntt4]=1
                    else:
                        Dmat_br[cntt1,cntt2,cntt3,cc,cntt4]=0                
Dmat_b=np.sum(Dmat_br,4)/17. 
#b2=30
#for cc in range(12):
#    b1=b2-10; b2=b1+20
#    for cntt1 in range(Dmat_b.shape[0]):
#        for cntt2 in range(Dmat_b.shape[1]):
#            for cntt3 in range(Dmat_b.shape[2]):                
#                Dmat_b[cntt1,cntt2,cntt3,cc]=np.sum(Dmat_br[cntt1,cntt2,cntt3,b1:b2,:]>0)/17.
print "prob_map fin"
import matplotlib.pyplot as plt

import copy
#l_new=[u'AnteriorTemporal', u'AnteriorTemporal', u'Ant.Fusiform', u'Ant.Fusiform', u'Ant.Sup.Frontal', u'Ant.Sup.Frontal', u'Bankssts', u'Bankssts', u'Caudal.Mid.Front.', u'Caudal.Mid.Front.', u'Cuneus', u'Cuneus', u'Frontalpole', u'Frontalpole', u'Inf.PostCentral', u'Inf.PostCentral', u'Inf.PreCentral', u'Inf.PreCentral', u'Inf.Parietal', u'Inf.Parietal', u'Inf.Temporal', u'Inf.Temporal', u'Lat.Occipital', u'Lat.Occipital', u'Lat.OrbitoFrontal', u'Lat.OrbitoFrontal', u'Lingual', u'Lingual', u'Mid.PostCentral', u'Mid.PostCentral', u'Mid.PreCentral', u'Mid.PreCentral', u'Mid.SuperiorFrontal', u'Mid.SuperiorFrontal', u'Mid.Temporal', u'Mid.Temporal', u'ParaCentral', u'ParaCentral', u'Parsopercularis', u'Parsopercularis', u'Parsorbitalis', u'Parsorbitalis', u'Parstriangularis', u'Parstriangularis', u'Pericalcarine', u'Pericalcarine', u'Post.Fusiform', u'Post.Fusiform', u'Post.SuperiorFrontal', u'Post.SuperiorFrontal', u'Precuneus', u'Precuneus', u'RostralMid.Frontal', u'RostralMid.Frontal', u'Sup.PostCentral', u'Sup.PostCentral', u'Sup.PreCentral', u'Sup.PreCentral', u'Sup.Parietal', u'Sup.Parietal', u'Sup.Temporal', u'Sup.Temporal', u'Supramarginal', u'Supramarginal']
l_new=[u'ATemporal', u'ATemporal', u'ASFront', u'ASFront', u'CaudMFront', u'CaudMFront', u'IPostCent', u'IPostCent', u'IPreCent', u'IPreCent', u'IParietal', u'IParietal', u'ITemporal', u'ITemporal', u'LatOcci', u'LatOcci', u'LatOrbFront', u'LatOrbFront', u'MPostCent', u'MPostCent', u'MPreCent', u'MPreCent', u'MSupFront', u'MSupFront', u'MTemporal', u'MTemporal', u'ParsOpercul', u'ParsOpercul', u'ParsOrbit', u'ParsOrbit', u'ParsTriangl', u'ParsTriangl', u'PSFrontal', u'PSFrontal', u'RostMFront', u'RostMFront', u'SPostCent', u'SPostCent', u'SPreCent', u'SPreCent',  u'STemporal', u'STemporal',u'SParietal', u'SParietal', u'SupraMarg', u'SupraMarg']

#cmap=mne.viz.mne_analyze_colormap(limits=[2,7,15], format='matplotlib')
print "plotting"
for ii,trange in zip(range(Dmat_b.shape[3]),['-100-0ms','-50-50ms', '0-100ms', '50-150ms','100-200ms','150-250ms', '200-300ms', '250-350ms','300-400ms','350-450ms', '400-500ms', '450-550ms']):#['50-150ms','150-250ms','250-350ms','350-450ms','450-550ms']):#time Dmatc.shape[2]):
    for jj,bands in zip(range(Dmat_b.shape[2]),['Theta','Alpha','Beta','Gamma']):#freq bands
        label_names1=copy.deepcopy(l_new)
        for kk in range(len(label_names)):
            posl=(Dmat_b[:,kk,jj,ii]>=16/17.).sum()
            if label_names[kk].endswith('lh'):
                label_names1[kk]=l_new[kk].replace(l_new[kk].encode(),l_new[kk].encode()+' ('+str(posl)+')')
            else:
                label_names1[kk]=l_new[kk].replace(l_new[kk].encode(),'('+str(posl)+') '+l_new[kk].encode())
                
        #dg=Dmat_b[Dmat_b[:,:,jj,ii].nonzero()[0],Dmat_b[:,:,jj,ii].nonzero()[1],jj,ii]; indices=Dmat_b[:,:,jj,ii].nonzero()
        #plot_connectivity_circle(dg,node_names=label_names1, indices=indices, node_angles=node_angles, node_colors=label_colors, vmin=-np.max(np.abs(dg)),vmax=np.max(np.abs(dg)),  colormap=cmap, colorbar=True,colorbar_pos=(-3.3,0.06),colorbar_size=0.15,fontsize_names=10,fontsize_colorbar=5,facecolor='black',textcolor='white')#,#######subplot=(5, 4, ii*4 + jj + 1),,####Dmatc[Dmatc[:,:,ii,2].nonzero()[0],Dmatc[:,:,ii,2].nonzero()[1],ii,2], indices=Dmatc[:,:,ii,2].nonzero(),padding=0, fontsize_colorbar=10, fig=fig,title='GrandAverage ppc Concrete-Baseline '+bands+' '+trange,
        plot_connectivity_circle(Dmat_b[:,:,jj,ii],node_names=label_names1,  node_angles=node_angles, node_colors=label_colors, vmin=0.3,vmax=1,  colormap='hot', colorbar=True,colorbar_pos=(-3.3,0.06),colorbar_size=0.15,fontsize_names=10,fontsize_colorbar=5,facecolor='black',textcolor='white')#,#######subplot=(5, 4, ii*4 + jj + 1),,####Dmatc[Dmatc[:,:,ii,2].nonzero()[0],Dmatc[:,:,ii,2].nonzero()[1],ii,2], indices=Dmatc[:,:,ii,2].nonzero(),padding=0, fontsize_colorbar=10, fig=fig,title='GrandAverage ppc Concrete-Baseline '+bands+' '+trange,

        out_name=out_path+'ppc_probmap_both_'+bands+'_'+trange+'.jpg'
        plt.savefig(out_name, facecolor='black',format='jpg',dpi=150)#, textcolor='white')#,tight_layout=True)
        plt.close()
        
out_mat=out_path+'ppc_Dmat_b.mat'
scipy.io.savemat(out_mat,{'Dmat_b':Dmat_b})