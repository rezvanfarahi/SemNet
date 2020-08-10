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
Dmat_cncrt=np.zeros((17,64,4,1201))
Dmat_abs=np.zeros((17,64,4,1201))
Dmat_sub=np.zeros((17,64,4,1201))

for cnt, meg in enumerate(ll):
    
    print cnt
    
    path_c=data_path+ meg+'degree_contrast_cb_absolute_pca_ppc.mat'
    ltc=scipy.io.loadmat(path_c,mat_dtype=True); deg_c=ltc['deg_c']
    Dmat_cncrt[cnt,:,:,:]=deg_c
    path_a=data_path+ meg+'degree_contrast_ab_absolute_pca_ppc.mat'
    lta=scipy.io.loadmat(path_a,mat_dtype=True); deg_a=lta['deg_a']
    Dmat_abs[cnt,:,:,:]=deg_a
    path_s=data_path+ meg+'degree_contrast_ca_absolute_pca_ppc.mat'
    lts=scipy.io.loadmat(path_s,mat_dtype=True); deg_s=lts['deg_s']
    Dmat_sub[cnt,:,:,:]=deg_s

print "deg files loaded, resampling"
    
b1=0; b2=20;
Dc=np.zeros((17,64,4,120))
Da=np.zeros((17,64,4,120))
Ds=np.zeros((17,64,4,120))
t_s2=[];
for ii in range(Dmat_cncrt.shape[3]/10):
    Dc[:,:,:,ii]=np.mean(Dmat_cncrt[:,:,:,b1:b2],3)
    Da[:,:,:,ii]=np.mean(Dmat_abs[:,:,:,b1:b2],3)
    Ds[:,:,:,ii]=np.mean(Dmat_sub[:,:,:,b1:b2],3)

    b1=b2-10;b2=b1+20
Dc1=Dc[:,:,0,50:110].squeeze().transpose([0,2,1])
Da1=Da[:,:,0,50:110].squeeze().transpose([0,2,1])
Ds1=Ds[:,:,0,50:110].squeeze().transpose([0,2,1])

#tc,cc,cpc,hc=mne.stats.permutation_cluster_1samp_test(Dc1, tail=0, max_step=0)
#ta,ca,cpa,ha=mne.stats.permutation_cluster_1samp_test(Da1, tail=0, max_step=0)
#ts,cs,cps,hs=mne.stats.permutation_cluster_1samp_test(Ds1, tail=0, max_step=0)

tc,pc,hc=mne.stats.permutation_t_test(Dc1.reshape(Dc1.shape[0],Dc1.shape[1]*Dc1.shape[2]), tail=0)
tc=tc.reshape(Dc1.shape[1],Dc1.shape[2])
pc=pc.reshape(Dc1.shape[1],Dc1.shape[2]); pc=-np.log10(pc)
print "cncrt fin"
ta,pa,ha=mne.stats.permutation_t_test(Da1.reshape(Da1.shape[0],Da1.shape[1]*Da1.shape[2]), tail=0)
ta=ta.reshape(Da1.shape[1],Da1.shape[2])
pa=pa.reshape(Da1.shape[1],Da1.shape[2]); pa=-np.log10(pa)
print "abs fin"
ts,ps,hs=mne.stats.permutation_t_test(Ds1.reshape(Ds1.shape[0],Ds1.shape[1]*Ds1.shape[2]), tail=0)
ts=ts.reshape(Ds1.shape[1],Ds1.shape[2])
ps=ps.reshape(Ds1.shape[1],Ds1.shape[2]); ps=-np.log10(ps)
print "sbtr fin"
#T_obs_plotc = np.nan * np.ones_like(tc)
#for c, p_val in zip(cc, cpc):
#    if p_val <= 0.05:
#        T_obs_plotc[c] = tc[c]
#T_obs_plotc=T_obs_plotc.transpose()   
#tc=tc.transpose()       
#vmax = np.max(np.abs(tc))
#vmin = -vmax
#plt.imshow(tc, cmap=plt.cm.gray,
#           extent=[0, 60, 1, 64],
#           aspect='auto', origin='lower', vmin=vmin, vmax=vmax)
#plt.imshow(T_obs_plotc, cmap=plt.cm.RdBu_r,
#           extent=[0, 60, 1, 64],
#           aspect='auto', origin='lower', vmin=vmin, vmax=vmax)
#plt.colorbar()
#plt.xlabel('time (ms)')
#plt.ylabel('Frequency (Hz)')
#plt.show()
#
#T_obs_plota = np.nan * np.ones_like(ta)
#for c, p_val in zip(ca, cpa):
#    if p_val <= 0.05:
#        T_obs_plota[c] = ta[c]
#T_obs_plota=T_obs_plota.transpose()   
#ta=ta.transpose()       
#vmax = np.max(np.abs(ta))
#vmin = -vmax
#plt.imshow(ta, cmap=plt.cm.gray,
#           extent=[0, 60, 1, 64],
#           aspect='auto', origin='lower', vmin=vmin, vmax=vmax)
#plt.imshow(T_obs_plota, cmap=plt.cm.RdBu_r,
#           extent=[0, 60, 1, 64],
#           aspect='auto', origin='lower', vmin=vmin, vmax=vmax)
#plt.colorbar()
#plt.xlabel('time (ms)')
#plt.ylabel('Frequency (Hz)')
#plt.show()
#
#T_obs_plots = np.nan * np.ones_like(ts)
#for c, p_val in zip(cs, cps):
#    if p_val <= 0.05:
#        T_obs_plots[c] = ts[c]
#T_obs_plots=T_obs_plots.transpose()   
#ts=ts.transpose()       
#vmax = np.max(np.abs(ts))
#vmin = -vmax
#plt.imshow(ts, cmap=plt.cm.gray,
#           extent=[0, 60, 1, 64],
#           aspect='auto', origin='lower', vmin=vmin, vmax=vmax)
#plt.imshow(T_obs_plots, cmap=plt.cm.RdBu_r,
#           extent=[0, 60, 1, 64],
#           aspect='auto', origin='lower', vmin=vmin, vmax=vmax)
#plt.colorbar()
#plt.xlabel('time (ms)')
#plt.ylabel('Frequency (Hz)')
#plt.show()

 
pc=pc.transpose()       
vmax = np.max(np.abs(pc))
vmin = -vmax

plt.imshow(pc, cmap=plt.cm.RdBu_r,
           extent=[0, 60, 1, 64],
           aspect='auto', origin='lower', vmin=vmin, vmax=vmax)
plt.colorbar()
plt.xlabel('time (ms)')
plt.ylabel('Frequency (Hz)')
plt.show()

pa=pa.transpose()       
vmax = np.max(np.abs(pa))
vmin = -vmax

plt.imshow(pa, cmap=plt.cm.RdBu_r,
           extent=[0, 60, 1, 64],
           aspect='auto', origin='lower', vmin=vmin, vmax=vmax)
plt.colorbar()
plt.xlabel('time (ms)')
plt.ylabel('Frequency (Hz)')
plt.show()

ps=ps.transpose()       
vmax = np.max(np.abs(ps))
vmin = -vmax

plt.imshow(ps, cmap=plt.cm.RdBu_r,
           extent=[0, 60, 1, 64],
           aspect='auto', origin='lower', vmin=vmin, vmax=vmax)
plt.colorbar()
plt.xlabel('time (ms)')
plt.ylabel('Frequency (Hz)')
plt.show()