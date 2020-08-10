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
sys.path.insert(1,'/imaging/local/software/mne_python/latest_v0.9')
sys.path.insert(1,'/imaging/local/software/mne_python/latest_v0.9/scikit-learn')

# for qsub
# (add ! here if needed) 
sys.path.append('/imaging/local/software/anaconda/latest/x86_64/bin/python')
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
import matplotlib.pyplot as plt


data_path = '/imaging/rf02/TypLexMEG/' # root directory for your MEG data
subjects_dir = '/imaging/rf02/TypLexMEG/'    # where your MRI subdirectories are

event_path = event_path = '/group/erp/data/olaf.hauk/MEG/TypLexM/data/MF22_new/'   # where event files are

label_path = '/imaging/rf02/TypLexMEG/fsaverage/label/'

inv_path = '/imaging/rf02/TypLexMEG/'
inv_fname = 'InverseOperator_fft_1_48_clean_ica_EMEG-inv.fif'

# get indices for subjects to be processed from command line input
# 
print sys.argv
subject_inds = [0]
for ss in sys.argv[1:]:
   subject_inds.append( int( ss ) )



#subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]#0,1,2,3,4,5,6,7,8]#,
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

graphpath='icaanalysis_results/jpg/time_resolved/ROI_connectivity/graph/'
outgraphpath='icaanalysis_results/stc/GrandAverage/time_resolved/graph/'

for cnt, meg in enumerate(ll):
    
    print cnt
    fname_cncrt = data_path + meg + 'firstMorphed_SemLoc_icaclean_Concrete_Source_Evoked_m500_700'
    stc_cncrt = mne.read_source_estimate(fname_cncrt)
    path_c=data_path+ graphpath+'coh_concrete_hubs.mat'
    ltc=scipy.io.loadmat(path_c,mat_dtype=True); deg_c=ltc['deg_final']
    path_a=data_path+ graphpath+'coh_abstract_hubs.mat'
    lta=scipy.io.loadmat(path_a,mat_dtype=True); deg_a=lta['deg_final']
   
    fsave_vertices = [np.arange(10242), np.arange(10242)]
    Xct=np.zeros((20484, deg_c.shape[2]))
    Xca=np.zeros((20484, deg_c.shape[2]))
    Xcb=np.zeros((20484, deg_c.shape[2]))
    Xcg=np.zeros((20484, deg_c.shape[2]))
    Xat=np.zeros((20484, deg_a.shape[2]))
    Xaa=np.zeros((20484, deg_a.shape[2]))
    Xab=np.zeros((20484, deg_a.shape[2]))
    Xag=np.zeros((20484, deg_a.shape[2]))
    print "deg files loaded"
    labelss = mne.read_labels_from_annot(subject='fsaverage', parc='rezvan_parc_mne_ctf_mn', subjects_dir=data_path)
    labelss.remove(labelss[-2])
    labelss.remove(labelss[-1])
    for ccc,label in enumerate(labelss):  
        print ccc
        bb=stc_cncrt.in_label(label)
        if ccc%2==0:
            for cnl in bb.lh_vertno:
                Xct[cnl,:]=deg_c[ccc,0,:]
                Xca[cnl,:]=deg_c[ccc,1,:]
                Xcb[cnl,:]=deg_c[ccc,2,:]
                Xcg[cnl,:]=deg_c[ccc,3,:]
                
                Xat[cnl,:]=deg_a[ccc,0,:]
                Xaa[cnl,:]=deg_a[ccc,1,:]
                Xab[cnl,:]=deg_a[ccc,2,:]
                Xag[cnl,:]=deg_a[ccc,3,:]
        else: 
            for cnr in bb.rh_vertno:
                Xct[cnr+10242,:]=deg_c[ccc,0,:]
                Xca[cnr+10242,:]=deg_c[ccc,1,:]
                Xcb[cnr+10242,:]=deg_c[ccc,2,:]
                Xcg[cnr+10242,:]=deg_c[ccc,3,:]
                
                Xat[cnr+10242,:]=deg_a[ccc,0,:]
                Xaa[cnr+10242,:]=deg_a[ccc,1,:]
                Xab[cnr+10242,:]=deg_a[ccc,2,:]
                Xag[cnr+10242,:]=deg_a[ccc,3,:]
            
       
    tmin1=50
    tstep1=100
    vertices_to = [np.arange(10242), np.arange(10242)]
    Xct_stc = mne.SourceEstimate(Xct, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
    Xca_stc = mne.SourceEstimate(Xca, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
    Xcb_stc = mne.SourceEstimate(Xcb, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
    Xcg_stc = mne.SourceEstimate(Xcg, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
    
    Xat_stc = mne.SourceEstimate(Xat, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
    Xaa_stc = mne.SourceEstimate(Xaa, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
    Xab_stc = mne.SourceEstimate(Xab, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
    Xag_stc = mne.SourceEstimate(Xag, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
    out_file_Xct=data_path + outgraphpath+'firstMorphed_SemLoc_icaclean_hubs_proportional_15_pca_coh_concrete_50_550_100ms_theta'
    Xct_stc.save(out_file_Xct)
    
    out_file_Xca=data_path + outgraphpath+'firstMorphed_SemLoc_icaclean_hubs_proportional_15_pca_coh_concrete_50_550_100ms_alpha'
    Xca_stc.save(out_file_Xca)
    
    out_file_Xcb=data_path + outgraphpath+'firstMorphed_SemLoc_icaclean_hubs_proportional_15_pca_coh_concrete_50_550_100ms_beta'
    Xcb_stc.save(out_file_Xcb)
    
    out_file_Xcg=data_path + outgraphpath+'firstMorphed_SemLoc_icaclean_hubs_proportional_15_pca_coh_concrete_50_550_100ms_gamma'
    Xcg_stc.save(out_file_Xcg)
    
    out_file_Xat=data_path + outgraphpath+'firstMorphed_SemLoc_icaclean_hubs_proportional_15_pca_coh_abstract_50_550_100ms_theta'
    Xat_stc.save(out_file_Xat)
    
    out_file_Xaa=data_path + outgraphpath+'firstMorphed_SemLoc_icaclean_hubs_proportional_15_pca_coh_abstract_50_550_100ms_alpha'
    Xaa_stc.save(out_file_Xaa)
    
    out_file_Xab=data_path + outgraphpath+'firstMorphed_SemLoc_icaclean_hubs_proportional_15_pca_coh_abstract_50_550_100ms_beta'
    Xab_stc.save(out_file_Xab)
    
    out_file_Xag=data_path + outgraphpath+'firstMorphed_SemLoc_icaclean_hubs_proportional_15_pca_coh_abstract_50_550_100ms_gamma'
    Xag_stc.save(out_file_Xag)

#b1=0; b2=20;
#Dc=np.zeros((17,64,4,120))
#Da=np.zeros((17,64,4,120))
#Ds=np.zeros((17,64,4,120))
#t_s2=[];
#for ii in range(Dmat_cncrt.shape[3]/10):
#    Dc[:,:,:,ii]=np.mean(Dmat_cncrt[:,:,:,b1:b2],3)
#    Da[:,:,:,ii]=np.mean(Dmat_abs[:,:,:,b1:b2],3)
#    Ds[:,:,:,ii]=np.mean(Dmat_sub[:,:,:,b1:b2],3)
#    b1=b2-10;b2=b1+20
#Dc1=Dc[:,:,0,50:110].squeeze().transpose([0,2,1])
#Da1=Da[:,:,0,50:110].squeeze().transpose([0,2,1])
#Ds1=Ds[:,:,0,50:110].squeeze().transpose([0,2,1])

