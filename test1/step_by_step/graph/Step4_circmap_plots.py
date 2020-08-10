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
#sys.path.insert(1,'/imaging/local/software/anaconda/1.9.1/x86_64/envs/py3k')
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
graphpath='icaanalysis_results/jpg/time_resolved/ROI_connectivity/graph/'

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
Dmat_cncrt=np.zeros((46,46,4,10,17))
Dmat_abs=np.zeros((46,46,4,10,17))
Dmat_both=np.zeros((46,46,4,10,17))

Dmat_c=np.zeros((46,46,4,5,17))
Dmat_a=np.zeros((46,46,4,5,17))
Dmat_b=np.zeros((46,46,4,5,17))
MTc=np.zeros((46,46,4,5))
MTa=np.zeros((46,46,4,5))
MTb=np.zeros((46,46,4,5))

for cnt1,meg in enumerate(ll):

#    print meg
#    #con_cncrt=np.divide(con_cncrt,con_blcn)
#    path_c=data_path+ meg+'plv_concrete_cwtmorlet.npy'
#    plv_concrete_cwtmorlet=np.load(path_c)
#    scipy.io.savemat(path_c[:-4],{'plv_concrete_cwtmorlet':plv_concrete_cwtmorlet})
#    print "con finished"
#    path_a=data_path+ meg+'plv_abstract_cwtmorlet.npy'
#    plv_abstract_cwtmorlet=np.load(path_a)
#    scipy.io.savemat(path_a[:-4],{'plv_abstract_cwtmorlet':plv_abstract_cwtmorlet})
#    print "abs finished"
#    path_s=data_path+ meg+'plv_subtract_cwtmorlet.npy'
#    plv_subtract_cwtmorlet=np.load(path_s)
#    scipy.io.savemat(path_s[:-4],{'plv_subtract_cwtmorlet':plv_subtract_cwtmorlet})
#    print "sub finished"
#    #con_abs=np.divide(con_abs,con_blab)

    #cnt1=cnt1+1
    
    print cnt1
    print meg
    fname_inv = inv_path + meg + 'InverseOperator_fft_1_48_clean_ica_EMEG-inv.fif'
    fname_raw = data_path + meg + '/semloc_ssstf_fft_1_48_clean_ica_raw.fif'
    fname_event = event_path + meg + 'semloc_raw_ssstnew.txt'

    # Load data
    inverse_operator = read_inverse_operator(fname_inv)
    srcin='/imaging/rf02/TypLexMEG/fsaverage/bem/fsaverage-ico-5-src.fif'
    src = mne.read_source_spaces(srcin)#inverse_operator['src']
    raw = Raw(fname_raw)
    
 
#    label_ts_concrete = mne.extract_label_time_course(stcs_concrete, labels[:-2], src, mode='pca_flip',return_generator=False)
#    label_ts_abstract = mne.extract_label_time_course(stcs_abstract, labels[:-2], src, mode='pca_flip',return_generator=False)
    path_c=data_path+ meg+'thresholded_15_pca_ppc_concrete_resampled_10win.mat'
    ltc=scipy.io.loadmat(path_c,mat_dtype=True); con_cncrt=ltc['thresh_c']; con_cncrt=np.abs(con_cncrt)
    #scipy.io.savemat(path_c,{'label_ts_concrete':label_ts_concrete})
    path_a=data_path+ meg+'thresholded_15_pca_ppc_abstract_resampled_10win.mat'
    #scipy.io.savemat(path_a,{'label_ts_abstract':label_ts_abstract})
    lta=scipy.io.loadmat(path_a,mat_dtype=True); con_abs=lta['thresh_a']; con_abs=np.abs(con_abs)
    path_b=data_path+ meg+'thresholded_15_pca_ppc_both_resampled_10win.mat'
    #scipy.io.savemat(path_a,{'label_ts_abstract':label_ts_abstract})
    ltb=scipy.io.loadmat(path_b,mat_dtype=True); con_both=ltb['thresh_b']; con_both=np.abs(con_both)
    
    print "lbl ts"

    
#    b=5    
#    for ccc in range(4):
#        Dmat_c[:,:,:,ccc,:]=Dmat_cncrt[:,:,:,b+ccc,:]-Dmat_cncrt[:,:,:,1,:]
#        Dmat_a[:,:,:,ccc,:]=Dmat_abs[:,:,:,b+ccc,:]-Dmat_abs[:,:,:,1,:]
#        Dmat_b[:,:,:,ccc,:]=Dmat_b[:,:,:,b+ccc,:]
    
    Dmat_cncrt[:,:,:,:,cnt1]=con_cncrt
    Dmat_abs[:,:,:,:,cnt1]=con_abs
    Dmat_both[:,:,:,:,cnt1]=con_both#(con_cncrt+con_abs)/2
path_bf=data_path+ graphpath+'ppc_degrees_15_resampled_10win_bothcon2_hubs_prob.mat'
ltc=scipy.io.loadmat(path_bf,mat_dtype=True); deg_finald=ltc['deg_final']
path_bf=data_path+ graphpath+'ppc_betweenness_15_resampled_10win_bothcon2_hubs_prob.mat'
ltc=scipy.io.loadmat(path_bf,mat_dtype=True); deg_finalb=ltc['deg_final']
labels = mne.read_labels_from_annot(subject='fsaverage', parc='rezvan_parc_mne_ctf_mn', subjects_dir=subjects_dir)
labels=labels[:-2]
label_colors = [label.color for label in labels]

path_cmap=data_path+ 'icaanalysis_results/rez_cmap.mat'
cmapp=scipy.io.loadmat(path_cmap,mat_dtype=True); cm=cmapp['rez_cmap']
import matplotlib as mpl
cmap=mpl.colors.ListedColormap(cm/255)
# Now, we visualize the connectivity using a circular graph layout

# First, we reorder the labels based on their location in the left hemi
labels[0].name=labels[0].name.replace('ATLnew','anterior_temporal')
labels[1].name=labels[1].name.replace('ATLnew','anterior_temporal')
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
    
#b2=10
#for cc in range(5):
#    b1=b2; b2=b1+20
#    Dmat_c[:,:,:,cc,:]=Dmat_cncrt[:,:,:,b1:b2,:].mean(axis=3)#-Dmat_cncrt[:,:,:,80:100,:].mean(axis=3)
#    Dmat_a[:,:,:,cc,:]=Dmat_abs[:,:,:,b1:b2,:].mean(axis=3)#-Dmat_abs[:,:,:,80:100,:].mean(axis=3)
#    Dmat_b[:,:,:,cc,:]=Dmat_both[:,:,:,b1:b2,:].mean(axis=3)
Dmat_c=Dmat_cncrt.copy()
print Dmat_c.shape
Dmat_a=Dmat_abs.copy()
Dmat_b=Dmat_both.copy()
Dmat_b[Dmat_b>0]=1
    
#for cntt1 in range(Dmat_c.shape[0]):
#    for cntt2 in range(Dmat_c.shape[1]):
#        if cntt2>cntt1:
#            Dmat_c[cntt1,cntt2,:,:,:]=Dmat_c[cntt2,cntt1,:,:,:]
#            Dmat_a[cntt1,cntt2,:,:,:]=Dmat_a[cntt2,cntt1,:,:,:]
#            Dmat_b[cntt1,cntt2,:,:,:]=Dmat_b[cntt2,cntt1,:,:,:]
#    

#Dmat_bz=np.zeros(Dmat_b.shape)
#for jj in range(Dmat_b.shape[2]): #bands
#    print "jj"
#    for ii in range(Dmat_b.shape[3]):#time
#        for kk in range(Dmat_b.shape[4]):#subjs
#            windeg=Dmat_b[:,:,jj,ii,kk];
#            Dmat_bz[:,:,jj,ii,kk]=(Dmat_b[:,:,jj,ii,kk]-windeg.mean())/windeg.std();
#Dmat_bzb=Dmat_bz.copy()
#Dmat_bzb[Dmat_bzb<1]=0
#Dmatc=Dmat_c.mean(axis=4)
#Dmata=Dmat_a.mean(axis=4)
#Dmatb=Dmat_bz.mean(axis=4)
####Dmatb[Dmatb<1]=0
#dd=Dmatc[Dmatc.nonzero()]
#thresh=sorted(np.abs(dd),reverse=True)[int(len(dd)/5)]
#Dmatc[np.abs(Dmatc)<thresh]=0
#
#dd=Dmata[Dmata.nonzero()]
#thresh=sorted(np.abs(dd),reverse=True)[int(len(dd)/5)]
#Dmata[np.abs(Dmata)<thresh]=0

#dd=Dmatb[Dmatb.nonzero()]
#thresh=sorted(dd,reverse=True)[int(len(dd)/5)]
#Dmatb[np.abs(Dmatb)<thresh]=0

#Dmat_bzp=np.zeros(Dmat_b.shape[:-1])
#for ii in range(Dmat_bzp.shape[3]): #time
#    for nn1 in range(Dmat_bzp.shape[0]):#nodes1
#        for nn2 in range(Dmat_bzp.shape[1]):
#            for bnd in range(Dmat_bzp.shape[2]): #bands
#                windeg=np.squeeze(Dmat_b[nn1,nn2,bnd,ii,:]);
#                Dmat_bzp[nn1,nn2,bnd,ii]=np.sum(windeg)/len(windeg)
#Dmatb=Dmat_bzp.copy()
#Dmatb[Dmatb<Dmatb.mean()+Dmatb.std()]=0
Dmatb=Dmat_b.mean(axis=4)
Dmatb[Dmatb<0.1]=0
Dmatb=Dmatb[:,:,:,1::2]


	# Plot the graph using node colors from the FreeSurfer parcellation. We only
	# show the 300 strongest connections.
#plot_connectivity_circle(con_res[:,:,0], label_names, n_lines=300,node_angles=node_angles, node_colors=label_colors,title='All-to-All Connectivity Concrete Word' 'Condition (ppc)')
print "averaging fin"
import matplotlib.pyplot as plt
	#plt.savefig('circle.png', facecolor='black')

	# Plot connectivity for both methods in the same plot
#no_names = [''] * len(label_names)

import copy
#l_new=[u'AnteriorTemporal', u'AnteriorTemporal', u'Ant.Fusiform', u'Ant.Fusiform', u'Ant.Sup.Frontal', u'Ant.Sup.Frontal', u'Bankssts', u'Bankssts', u'Caudal.Mid.Front.', u'Caudal.Mid.Front.', u'Cuneus', u'Cuneus', u'Frontalpole', u'Frontalpole', u'Inf.PostCentral', u'Inf.PostCentral', u'Inf.PreCentral', u'Inf.PreCentral', u'Inf.Parietal', u'Inf.Parietal', u'Inf.Temporal', u'Inf.Temporal', u'Lat.Occipital', u'Lat.Occipital', u'Lat.OrbitoFrontal', u'Lat.OrbitoFrontal', u'Lingual', u'Lingual', u'Mid.PostCentral', u'Mid.PostCentral', u'Mid.PreCentral', u'Mid.PreCentral', u'Mid.SuperiorFrontal', u'Mid.SuperiorFrontal', u'Mid.Temporal', u'Mid.Temporal', u'ParaCentral', u'ParaCentral', u'Parsopercularis', u'Parsopercularis', u'Parsorbitalis', u'Parsorbitalis', u'Parstriangularis', u'Parstriangularis', u'Pericalcarine', u'Pericalcarine', u'Post.Fusiform', u'Post.Fusiform', u'Post.SuperiorFrontal', u'Post.SuperiorFrontal', u'Precuneus', u'Precuneus', u'RostralMid.Frontal', u'RostralMid.Frontal', u'Sup.PostCentral', u'Sup.PostCentral', u'Sup.PreCentral', u'Sup.PreCentral', u'Sup.Parietal', u'Sup.Parietal', u'Sup.Temporal', u'Sup.Temporal', u'Supramarginal', u'Supramarginal']
#l_new=[u'ATemporal', u'ATemporal', u'AIFront', u'AIFront', u'ASFront', u'ASFront', u'CaudMFront', u'CaudMFront', u'IPostCent', u'IPostCent', u'IPreCent', u'IPreCent', u'IParietal', u'IParietal', u'ITemporal', u'ITemporal', u'LatOcci', u'LatOcci', u'MPostCent', u'MPostCent', u'MPreCent', u'MPreCent', u'MSupFront', u'MSupFront', u'MTemporal', u'MTemporal', u'ParaCent', u'ParaCent', u'PFusiform', u'PFusiform', u'PIFrontal', u'PIFrontal',u'PSFrontal', u'PSFrontal', u'RostMFront', u'RostMFront', u'SPostCent', u'SPostCent', u'SPreCent', u'SPreCent',  u'STemporal', u'STemporal',u'SParietal', u'SParietal', u'SupraMarg', u'SupraMarg']
l_new=[u'ATL', u'ATL', u'aIFG', u'aIFG', u'aSFG', u'aSFG', u'CMF', u'CMF', u'iPosC', u'iPosC', u'iPreC', u'iPreC', u'IPL', u'IPL', u'ITL', u'ITL', u'LOcci', u'LOcci', u'mPosC', u'mPosC', u'mPreC', u'mPreC', u'mSFG', u'mSFG', u'MTL', u'MTL', u'ParaC', u'ParaC', u'pFusi', u'pFusi', u'pIFG', u'pIFG',u'pSFG', u'pSFG', u'RMF', u'RMF', u'SPosC', u'SPosC', u'SPreC', u'SPreC',  u'STL', u'STL',u'SPL', u'SPL', u'SMarg', u'SMarg']


#cmap=mne.viz.mne_analyze_colormap(limits=[2,7,15], format='matplotlib')
Dmatb1=Dmatb.copy()
for ii in range(deg_finald.shape[0]):
    for jj in range(deg_finald.shape[1]):
        for kk in range(deg_finald.shape[2]):
            for cc in range(Dmatb.shape[0]):
                if (Dmatb1[ii,cc,jj,kk]>0 and deg_finald[ii,jj,kk]==0):
                    Dmatb[ii,cc,jj,kk]=0.2#0.2
                    #Dmatb[cc,ii,jj,kk]=0.2#0.2
                if (Dmatb1[ii,cc,jj,kk]>0 and deg_finald[ii,jj,kk]>0):
                    Dmatb[ii,cc,jj,kk]=0.8#0.5
                    #Dmatb[cc,ii,jj,kk]=0.8#0.5
                
#                if (Dmatb1[ii,cc,jj,kk]>0 and deg_finald[ii,jj,kk]==0 and deg_finalb[ii,jj,kk]==0):
#                    Dmatb[ii,cc,jj,kk]=-0.4#0.2
#                    Dmatb[cc,ii,jj,kk]=-0.4#0.2
#                if (Dmatb1[ii,cc,jj,kk]>0 and np.logical_xor(deg_finald[ii,jj,kk]>0,deg_finalb[ii,jj,kk]>0)):
#                    Dmatb[ii,cc,jj,kk]=0.5#0.5
#                    Dmatb[cc,ii,jj,kk]=0.5#0.5
#                if (Dmatb1[ii,cc,jj,kk]>0 and (deg_finald[ii,jj,kk]>0 and deg_finalb[ii,jj,kk]>0)):
#                    Dmatb[ii,cc,jj,kk]=1.0#0.8
#                    Dmatb[cc,ii,jj,kk]=1.0#0.8
                
#                else:
#                    Dmatb[ii,cc,jj,kk]=0
#                    Dmatb[cc,ii,jj,kk]=0

                
        
for ii,trange in zip(range(Dmatb.shape[3]),['50-150ms','150-250ms', '250-350ms', '350-450ms', '450-550ms']):#['50-150ms','150-250ms','250-350ms','350-450ms','450-550ms']):#time Dmatc.shape[2]):
    for jj,bands in zip(range(Dmatb.shape[2]),['Theta','Alpha','Beta','Gamma']):#freq bands
        label_names1=copy.deepcopy(l_new)
        for kk in range(len(label_names)):
            if deg_finald[kk,jj,ii]>0:
                posl='*'
            else:
                posl=''
            #posl=(Dmatb[:,kk,jj,ii]==1).sum()
#            if np.logical_xor(deg_finald[kk,jj,ii]>0,deg_finalb[kk,jj,ii]>0):
#                posl='*'
#            if (deg_finald[kk,jj,ii]>0 and deg_finalb[kk,jj,ii]>0):
#                posl='**'
#            if (deg_finald[kk,jj,ii]==0 and deg_finalb[kk,jj,ii]==0):
#                posl=''
#            if label_names[kk].endswith('lh'):
#                label_names1[kk]=l_new[kk].replace(l_new[kk].encode(),l_new[kk].encode()+' ('+str(posl)+')')
#            else:
#                label_names1[kk]=l_new[kk].replace(l_new[kk].encode(),'('+str(posl)+') '+l_new[kk].encode())
#              
            if label_names[kk].endswith('lh'):
                label_names1[kk]=l_new[kk].replace(l_new[kk].encode(),posl+l_new[kk].encode()+posl)
            else:
                label_names1[kk]=l_new[kk].replace(l_new[kk].encode(),posl+l_new[kk].encode()+posl)
#            else:
#                if label_names[kk].endswith('lh'):
#                    label_names1[kk]=l_new[kk].replace(l_new[kk].encode(),l_new[kk].encode())
#                else:
#                    label_names1[kk]=l_new[kk].replace(l_new[kk].encode(),l_new[kk].encode())
              
              
        #dg=Dmatb[Dmatb[:,:,jj,ii].nonzero()[0],Dmatb[:,:,jj,ii].nonzero()[1],jj,ii]; indices=Dmatb[:,:,jj,ii].nonzero()
        plot_connectivity_circle(Dmatb[:,:,jj,ii],node_names=label_names1, node_angles=node_angles, node_colors=label_colors, vmin=0,vmax=1,  colormap='hot', colorbar=False,colorbar_pos=(-3.3,0.06),colorbar_size=0.15,fontsize_names=14,fontsize_colorbar=5,facecolor='black',textcolor='white')#facecolor='white',textcolor='black',#######subplot=(5, 4, ii*4 + jj + 1),,####Dmatc[Dmatc[:,:,ii,2].nonzero()[0],Dmatc[:,:,ii,2].nonzero()[1],ii,2], indices=Dmatc[:,:,ii,2].nonzero(),padding=0, fontsize_colorbar=10, fig=fig,title='GrandAverage ppc Concrete-Abstract '+bands+' '+trange, 
        out_name=out_path+'ppc_15_resampled_10win_bothcon_hubs_prob_'+bands+'_'+trange+'.jpg'
        plt.savefig(out_name, facecolor='black',format='jpg',dpi=150)#,tight_layout=True)
        plt.close()
        
#plt.show()

