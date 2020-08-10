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
out_path='/imaging/rf02/TypLexMEG/icaanalysis_results/jpg/time_resolved/ROI_connectivity/new/ppc/'

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
Dmat_cncrt=np.zeros((64,64,4,1201,17))
Dmat_abs=np.zeros((64,64,4,1201,17))
Dmat_sub=np.zeros((64,64,4,1201,17))

Dmat_c=np.zeros((64,64,4,4,17))
Dmat_a=np.zeros((64,64,4,4,17))
Dmat_s=np.zeros((64,64,4,4,17))
MTc=np.zeros((64,64,4,4))
MTa=np.zeros((64,64,4,4))
MTs=np.zeros((64,64,4,4))

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
    path_c=data_path+ meg+'pca_ppc_concrete_cwtmorlet.mat'
    ltc=scipy.io.loadmat(path_c,mat_dtype=True); con_cncrt=ltc['con_cncrt']; con_cncrt=np.abs(con_cncrt)
    #scipy.io.savemat(path_c,{'label_ts_concrete':label_ts_concrete})
    path_a=data_path+ meg+'pca_ppc_abstract_cwtmorlet.mat'
    #scipy.io.savemat(path_a,{'label_ts_abstract':label_ts_abstract})
    lta=scipy.io.loadmat(path_a,mat_dtype=True); con_abs=lta['con_abs']; con_abs=np.abs(con_abs)
    print "lbl ts"

    labels = mne.read_labels_from_annot(subject='fsaverage', parc='rezvan_annot', subjects_dir=subjects_dir)
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
#    b=5    
#    for ccc in range(4):
#        Dmat_c[:,:,:,ccc,:]=Dmat_cncrt[:,:,:,b+ccc,:]-Dmat_cncrt[:,:,:,1,:]
#        Dmat_a[:,:,:,ccc,:]=Dmat_abs[:,:,:,b+ccc,:]-Dmat_abs[:,:,:,1,:]
#        Dmat_s[:,:,:,ccc,:]=Dmat_s[:,:,:,b+ccc,:]
    
    Dmat_cncrt[:,:,:,:,cnt1]=con_cncrt
    Dmat_abs[:,:,:,:,cnt1]=con_abs
    Dmat_sub[:,:,:,:,cnt1]=con_cncrt-con_abs
    
b2=649
for cc in range(4):
    b1=b2-100+1; b2=b1+200
    Dmat_c[:,:,:,cc,:]=Dmat_cncrt[:,:,:,b1:b2,:].mean(axis=3)-Dmat_cncrt[:,:,:,200:400,:].mean(axis=3)
    Dmat_a[:,:,:,cc,:]=Dmat_abs[:,:,:,b1:b2,:].mean(axis=3)-Dmat_abs[:,:,:,200:400,:].mean(axis=3)
    Dmat_s[:,:,:,cc,:]=Dmat_sub[:,:,:,b1:b2,:].mean(axis=3)
    
for cntt1 in range(Dmat_c.shape[0]):
    for cntt2 in range(Dmat_c.shape[1]):
        if cntt2>cntt1:
            Dmat_c[cntt1,cntt2]=Dmat_c[cntt2,cntt1]
            Dmat_a[cntt1,cntt2]=Dmat_a[cntt2,cntt1]
            Dmat_s[cntt1,cntt2]=Dmat_s[cntt2,cntt1]
    

Dmatc=Dmat_c.mean(axis=4)
Dmata=Dmat_a.mean(axis=4)
Dmats=Dmat_s.mean(axis=4)

dd=Dmatc[Dmatc.nonzero()]
thresh=sorted(np.abs(dd),reverse=True)[int(len(dd)/5)]
Dmatc[np.abs(Dmatc)<thresh]=0

dd=Dmata[Dmata.nonzero()]
thresh=sorted(np.abs(dd),reverse=True)[int(len(dd)/5)]
Dmata[np.abs(Dmata)<thresh]=0

dd=Dmats[Dmats.nonzero()]
thresh=sorted(np.abs(dd),reverse=True)[int(len(dd)/5)]
Dmats[np.abs(Dmats)<thresh]=0
	# Plot the graph using node colors from the FreeSurfer parcellation. We only
	# show the 300 strongest connections.
#plot_connectivity_circle(con_res[:,:,0], label_names, n_lines=300,node_angles=node_angles, node_colors=label_colors,title='All-to-All Connectivity Concrete Word' 'Condition (ppc)')
print "averaging fin"
import matplotlib.pyplot as plt
	#plt.savefig('circle.png', facecolor='black')

	# Plot connectivity for both methods in the same plot
#no_names = [''] * len(label_names)
for cc1 in range(Dmat_c.shape[3]):#time
    for cc2 in range(Dmat_c.shape[2]):#band
        for cc3 in range(Dmat_c.shape[1]):#conn
            for cc4 in range(Dmat_c.shape[0]):#conn
                MTc[cc4,cc3,cc2,cc1]= mne.stats.ttest_1samp_no_p(Dmat_c[cc4,cc3,cc2,cc1,:], sigma=1e-3, method='relative');
                MTa[cc4,cc3,cc2,cc1]= mne.stats.ttest_1samp_no_p(Dmat_a[cc4,cc3,cc2,cc1,:], sigma=1e-3, method='relative');
                MTs[cc4,cc3,cc2,cc1]= mne.stats.ttest_1samp_no_p(Dmat_s[cc4,cc3,cc2,cc1,:], sigma=1e-3, method='relative');
MTc[np.isnan(MTc)]=0
MTa[np.isnan(MTa)]=0
MTs[np.isnan(MTs)]=0
          
print "ttest fin"
            
#Matxc=Dmat_c.reshape(Dmat_c.shape[4],Dmat_c.shape[0]*Dmat_c.shape[1]*Dmat_c.shape[2]*Dmat_c.shape[3])
#Matxa=Dmat_a.reshape(Dmat_a.shape[4],Dmat_a.shape[0]*Dmat_a.shape[1]*Dmat_a.shape[2]*Dmat_a.shape[3])
#Matxs=Dmat_s.reshape(Dmat_s.shape[4],Dmat_s.shape[0]*Dmat_s.shape[1]*Dmat_s.shape[2]*Dmat_s.shape[3])
p_threshold = 0.05
n_subjects=17
t_threshold = -scipy.stats.distributions.t.ppf(p_threshold/2., n_subjects - 1)
#MTc= mne.stats.ttest_1samp_no_p(Matxc, sigma=0, method='relative');MTc=MTc.reshape(Dmatc.shape)
MTc[np.abs(MTc)<np.abs(t_threshold)]=0
#MTa= mne.stats.ttest_1samp_no_p(Matxa, sigma=0, method='relative');MTa=MTa.reshape(Dmata.shape)
MTa[np.abs(MTa)<np.abs(t_threshold)]=0
#MTs= mne.stats.ttest_1samp_no_p(Matxs, sigma=0, method='relative');MTs=MTs.reshape(Dmats.shape)
MTs[np.abs(MTs)<np.abs(t_threshold)]=0
import copy
#l_new=[u'AnteriorTemporal', u'AnteriorTemporal', u'Ant.Fusiform', u'Ant.Fusiform', u'Ant.Sup.Frontal', u'Ant.Sup.Frontal', u'Bankssts', u'Bankssts', u'Caudal.Mid.Front.', u'Caudal.Mid.Front.', u'Cuneus', u'Cuneus', u'Frontalpole', u'Frontalpole', u'Inf.PostCentral', u'Inf.PostCentral', u'Inf.PreCentral', u'Inf.PreCentral', u'Inf.Parietal', u'Inf.Parietal', u'Inf.Temporal', u'Inf.Temporal', u'Lat.Occipital', u'Lat.Occipital', u'Lat.OrbitoFrontal', u'Lat.OrbitoFrontal', u'Lingual', u'Lingual', u'Mid.PostCentral', u'Mid.PostCentral', u'Mid.PreCentral', u'Mid.PreCentral', u'Mid.SuperiorFrontal', u'Mid.SuperiorFrontal', u'Mid.Temporal', u'Mid.Temporal', u'ParaCentral', u'ParaCentral', u'Parsopercularis', u'Parsopercularis', u'Parsorbitalis', u'Parsorbitalis', u'Parstriangularis', u'Parstriangularis', u'Pericalcarine', u'Pericalcarine', u'Post.Fusiform', u'Post.Fusiform', u'Post.SuperiorFrontal', u'Post.SuperiorFrontal', u'Precuneus', u'Precuneus', u'RostralMid.Frontal', u'RostralMid.Frontal', u'Sup.PostCentral', u'Sup.PostCentral', u'Sup.PreCentral', u'Sup.PreCentral', u'Sup.Parietal', u'Sup.Parietal', u'Sup.Temporal', u'Sup.Temporal', u'Supramarginal', u'Supramarginal']
l_new=[u'ATemporal', u'ATemporal', u'AFusiform', u'AFusiform', u'ASFront', u'ASFront', u'Bankssts', u'Bankssts', u'CaudMFront', u'CaudMFront', u'Cuneus', u'Cuneus', u'FrontPole', u'FrontPole', u'IPostCent', u'IPostCent', u'IPreCent', u'IPreCent', u'IParietal', u'IParietal', u'ITemporal', u'ITemporal', u'LatOcci', u'LatOcci', u'LatOrbFront', u'LatOrbFront', u'Lingual', u'Lingual', u'MPostCent', u'MPostCent', u'MPreCent', u'MPreCent', u'MSupFront', u'MSupFront', u'MTemporal', u'MTemporal', u'ParaCent', u'ParaCent', u'ParsOpercul', u'ParsOpercul', u'ParsOrbit', u'ParsOrbit', u'ParsTriangl', u'ParsTriangl', u'PeriCalc', u'PeriCalc', u'PFusiform', u'PFusiform', u'PSFrontal', u'PSFrontal', u'Precuneus', u'Precuneus', u'RostMFront', u'RostMFront', u'SPostCent', u'SPostCent', u'SPreCent', u'SPreCent', u'SParietal', u'SParietal', u'STemporal', u'STemporal', u'SupraMarg', u'SupraMarg']

#cmap=mne.viz.mne_analyze_colormap(limits=[2,7,15], format='matplotlib')
print "plotting"
for ii,trange in zip(range(Dmatc.shape[3]),['50-250ms','150-350ms', '250-450ms', '350-550ms']):#['50-150ms','150-250ms','250-350ms','350-450ms','450-550ms']):#time Dmatc.shape[2]):
    for jj,bands in zip(range(Dmatc.shape[2]),['Theta','Alpha','Beta','Gamma']):#freq bands
        label_names1=copy.deepcopy(l_new)
        for kk in range(len(label_names)):
            posl=(Dmatc[:,kk,jj,ii]>0).sum()
            negl=(Dmatc[:,kk,jj,ii]<0).sum()
            if label_names[kk].endswith('lh'):
                label_names1[kk]=l_new[kk].replace(l_new[kk].encode(),l_new[kk].encode()+' (+'+str(posl)+', -'+str(negl)+')')
            else:
                label_names1[kk]=l_new[kk].replace(l_new[kk].encode(),'(+'+str(posl)+', -'+str(negl)+') '+l_new[kk].encode())
                
        dg=Dmatc[Dmatc[:,:,jj,ii].nonzero()[0],Dmatc[:,:,jj,ii].nonzero()[1],jj,ii]; indices=Dmatc[:,:,jj,ii].nonzero()
        plot_connectivity_circle(dg,node_names=label_names1, indices=indices, node_angles=node_angles, node_colors=label_colors, vmin=-np.max(np.abs(dg)),vmax=np.max(np.abs(dg)),  colormap=cmap, colorbar=True,colorbar_pos=(-3.3,0.06),colorbar_size=0.15,fontsize_names=10,fontsize_colorbar=5,facecolor='black',textcolor='white')#,#######subplot=(5, 4, ii*4 + jj + 1),,####Dmatc[Dmatc[:,:,ii,2].nonzero()[0],Dmatc[:,:,ii,2].nonzero()[1],ii,2], indices=Dmatc[:,:,ii,2].nonzero(),padding=0, fontsize_colorbar=10, fig=fig,title='GrandAverage ppc Concrete-Baseline '+bands+' '+trange,
        out_name=out_path+'ppc_grandaverage_cb_'+bands+'_'+trange+'.jpg'
        plt.savefig(out_name, facecolor='black',format='jpg',dpi=150)#, textcolor='white')#,tight_layout=True)
        plt.close()
        label_names1=copy.deepcopy(l_new)
        for kk in range(len(label_names)):
            posl=(MTc[:,kk,jj,ii]>0).sum()
            negl=(MTc[:,kk,jj,ii]<0).sum()
            if label_names[kk].endswith('lh'):
                label_names1[kk]=l_new[kk].replace(l_new[kk].encode(),l_new[kk].encode()+' (+'+str(posl)+', -'+str(negl)+')')
            else:
                label_names1[kk]=l_new[kk].replace(l_new[kk].encode(),'(+'+str(posl)+', -'+str(negl)+') '+l_new[kk].encode())
                
        dm=MTc[MTc[:,:,jj,ii].nonzero()[0],MTc[:,:,jj,ii].nonzero()[1],jj,ii]; indices=MTc[:,:,jj,ii].nonzero()
        plot_connectivity_circle(dm,node_names=label_names1, indices=indices, n_lines=None,node_angles=node_angles, node_colors=label_colors, vmin=-4, vmax=4, colormap=cmap, colorbar=True,colorbar_pos=(-3.3,0.06),colorbar_size=0.15,fontsize_names=10,fontsize_colorbar=6,facecolor='black',textcolor='white')#subplot=(5, 4, ii*4 + jj + 1),,####Dmatc[Dmatc[:,:,ii,2].nonzero()[0],Dmatc[:,:,ii,2].nonzero()[1],ii,2], indices=Dmatc[:,:,ii,2].nonzero(),padding=0, fontsize_colorbar=10, fig=fig,title='Univariate t-map ppc Concrete-Baseline '+bands+' '+trange,padding=10.0,colorbar_pos=(-0.45,0.14),
        out_name=out_path+'ppc_uvttest_cb_'+bands+'_'+trange+'.jpg'
        plt.savefig(out_name, facecolor='black',format='jpg',dpi=150)#,tight_layout=True)
        plt.close()
        
        
for ii,trange in zip(range(Dmata.shape[3]),['50-250ms','150-350ms', '250-450ms', '350-550ms']):#['50-150ms','150-250ms','250-350ms','350-450ms','450-550ms']):#time Dmatc.shape[2]):
    for jj,bands in zip(range(Dmata.shape[2]),['Theta','Alpha','Beta','Gamma']):#freq bands
        label_names1=copy.deepcopy(l_new)
        for kk in range(len(label_names)):
            posl=(Dmata[:,kk,jj,ii]>0).sum()
            negl=(Dmata[:,kk,jj,ii]<0).sum()
            if label_names[kk].endswith('lh'):
                label_names1[kk]=l_new[kk].replace(l_new[kk].encode(),l_new[kk].encode()+' (+'+str(posl)+', -'+str(negl)+')')
            else:
                label_names1[kk]=l_new[kk].replace(l_new[kk].encode(),'(+'+str(posl)+', -'+str(negl)+') '+l_new[kk].encode())
              
        dg=Dmata[Dmata[:,:,jj,ii].nonzero()[0],Dmata[:,:,jj,ii].nonzero()[1],jj,ii]; indices=Dmata[:,:,jj,ii].nonzero()
        plot_connectivity_circle(dg,node_names=label_names1, indices=indices, node_angles=node_angles, node_colors=label_colors, vmin=-np.max(np.abs(dg)),vmax=np.max(np.abs(dg)), colormap=cmap, colorbar=True,colorbar_pos=(-3.3,0.06),colorbar_size=0.15,fontsize_names=10,fontsize_colorbar=5,facecolor='black',textcolor='white')#facecolor='white',textcolor='black',#######subplot=(5, 4, ii*4 + jj + 1),,####Dmatc[Dmatc[:,:,ii,2].nonzero()[0],Dmatc[:,:,ii,2].nonzero()[1],ii,2], indices=Dmatc[:,:,ii,2].nonzero(),padding=0, fontsize_colorbar=10, fig=fig,title='GrandAverage ppc Abstract-Baseline '+bands+' '+trange,  
        out_name=out_path+'ppc_grandaverage_ab_'+bands+'_'+trange+'.jpg'
        plt.savefig(out_name, facecolor='black',format='jpg',dpi=150)#,tight_layout=True)
        plt.close()
        label_names1=copy.deepcopy(l_new)
        for kk in range(len(label_names)):
            posl=(MTa[:,kk,jj,ii]>0).sum()
            negl=(MTa[:,kk,jj,ii]<0).sum()
            if label_names[kk].endswith('lh'):
                label_names1[kk]=l_new[kk].replace(l_new[kk].encode(),l_new[kk].encode()+' (+'+str(posl)+', -'+str(negl)+')')
            else:
                label_names1[kk]=l_new[kk].replace(l_new[kk].encode(),'(+'+str(posl)+', -'+str(negl)+') '+l_new[kk].encode())
                
        dm=MTa[MTa[:,:,jj,ii].nonzero()[0],MTa[:,:,jj,ii].nonzero()[1],jj,ii]; indices=MTa[:,:,jj,ii].nonzero()
        plot_connectivity_circle(dm,node_names=label_names1, indices=indices, n_lines=None,node_angles=node_angles, node_colors=label_colors, vmin=-4, vmax=4, colormap=cmap, colorbar=True,colorbar_pos=(-3.3,0.06),colorbar_size=0.15,fontsize_names=10,fontsize_colorbar=6,facecolor='black',textcolor='white')#subplot=(5, 4, ii*4 + jj + 1),,####Dmatc[Dmatc[:,:,ii,2].nonzero()[0],Dmatc[:,:,ii,2].nonzero()[1],ii,2], indices=Dmatc[:,:,ii,2].nonzero(),padding=0, fontsize_colorbar=10, fig=fig,title='Univariate t-map ppc Abstract-Baseline '+bands+' '+trange,  
        out_name=out_path+'ppc_uvttest_ab_'+bands+'_'+trange+'.jpg'
        plt.savefig(out_name, facecolor='black',format='jpg',dpi=150)#,tight_layout=True)
        plt.close()
#Dmats=Dmatc-Dmata;
        
for ii,trange in zip(range(Dmats.shape[3]),['50-250ms','150-350ms', '250-450ms', '350-550ms']):#['50-150ms','150-250ms','250-350ms','350-450ms','450-550ms']):#time Dmatc.shape[2]):
    for jj,bands in zip(range(Dmats.shape[2]),['Theta','Alpha','Beta','Gamma']):#freq bands
        label_names1=copy.deepcopy(l_new)
        for kk in range(len(label_names)):
            posl=(Dmats[:,kk,jj,ii]>0).sum()
            negl=(Dmats[:,kk,jj,ii]<0).sum()
            if label_names[kk].endswith('lh'):
                label_names1[kk]=l_new[kk].replace(l_new[kk].encode(),l_new[kk].encode()+' (+'+str(posl)+', -'+str(negl)+')')
            else:
                label_names1[kk]=l_new[kk].replace(l_new[kk].encode(),'(+'+str(posl)+', -'+str(negl)+') '+l_new[kk].encode())
              
        dg=Dmats[Dmats[:,:,jj,ii].nonzero()[0],Dmats[:,:,jj,ii].nonzero()[1],jj,ii]; indices=Dmats[:,:,jj,ii].nonzero()
        plot_connectivity_circle(dg,node_names=label_names1, indices=indices, node_angles=node_angles, node_colors=label_colors, vmin=-np.max(np.abs(dg)),vmax=np.max(np.abs(dg)),  colormap=cmap, colorbar=True,colorbar_pos=(-3.3,0.06),colorbar_size=0.15,fontsize_names=10,fontsize_colorbar=5,facecolor='black',textcolor='white')#facecolor='white',textcolor='black',#######subplot=(5, 4, ii*4 + jj + 1),,####Dmatc[Dmatc[:,:,ii,2].nonzero()[0],Dmatc[:,:,ii,2].nonzero()[1],ii,2], indices=Dmatc[:,:,ii,2].nonzero(),padding=0, fontsize_colorbar=10, fig=fig,title='GrandAverage ppc Concrete-Abstract '+bands+' '+trange, 
        out_name=out_path+'ppc_grandaverage_ca_'+bands+'_'+trange+'.jpg'
        plt.savefig(out_name, facecolor='black',format='jpg',dpi=150)#,tight_layout=True)
        plt.close()
        label_names1=copy.deepcopy(l_new)
        for kk in range(len(label_names)):
            posl=(MTs[:,kk,jj,ii]>0).sum()
            negl=(MTs[:,kk,jj,ii]<0).sum()
            if label_names[kk].endswith('lh'):
                label_names1[kk]=l_new[kk].replace(l_new[kk].encode(),l_new[kk].encode()+' (+'+str(posl)+', -'+str(negl)+')')
            else:
                label_names1[kk]=l_new[kk].replace(l_new[kk].encode(),'(+'+str(posl)+', -'+str(negl)+') '+l_new[kk].encode())
                
        dm=MTs[MTs[:,:,jj,ii].nonzero()[0],MTs[:,:,jj,ii].nonzero()[1],jj,ii]; indices=MTs[:,:,jj,ii].nonzero()
        plot_connectivity_circle(dm,node_names=label_names1, indices=indices, n_lines=None,node_angles=node_angles, node_colors=label_colors, vmin=-4, vmax=4, colormap=cmap, colorbar=True,colorbar_pos=(-3.3,0.06),colorbar_size=0.15,fontsize_names=10,fontsize_colorbar=6,facecolor='black',textcolor='white')#subplot=(5, 4, ii*4 + jj + 1),,####Dmatc[Dmatc[:,:,ii,2].nonzero()[0],Dmatc[:,:,ii,2].nonzero()[1],ii,2], indices=Dmatc[:,:,ii,2].nonzero(),padding=0, fontsize_colorbar=10, fig=fig,title='Univariate t-map ppc Concrete-Abstract '+bands+' '+trange, 
        out_name=out_path+'ppc_uvttest_ca_'+bands+'_'+trange+'.jpg'
        plt.savefig(out_name, facecolor='black',format='jpg',dpi=150)#,tight_layout=True)
        plt.close()
#plt.show()

