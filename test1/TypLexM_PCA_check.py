# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 18:28:35 2015

@author: rf02
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 03:52:45 2015

@author: rf02
"""
print(__doc__)

import sys
sys.path.append('/imaging/local/software/python_packages/nibabel/1.3.0')
sys.path.append('/imaging/local/software/python_packages/pysurfer/v0.3.1')
sys.path.append('/imaging/local/software/python_packages/scikit-learn/v0.14.1')
sys.path.append('/imaging/local/software/mne_python/latest')

# for qsub
# (add ! here if needed) 
sys.path.append('/imaging/local/software/anaconda/latest/x86_64/bin/python')
sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
###

import numpy as np
import mne
import scipy


###############################################################################
out_path = '/imaging/rf02/TypLexMEG/icaanalysis_results/stc/GrandAverage/power/' # root directory for your MEG data
subjects_dir = '/imaging/rf02/TypLexMEG/'    # where your MRI subdirectories are
# where event-files are
data_path = '/imaging/rf02/TypLexMEG/'    # where event files are

label_path = '/home/rf02/rezvan/TypLexMEG/createdlabels_SL/'

inv_path = '/imaging/rf02/TypLexMEG/'
inv_fname = 'InverseOperator_EMEG-inv.fif'

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



# Labels/ROIs
#labellist = ['atlleft-lh']#, 'atlleftt-rh'
tmin, tmax = -0.5, 0.7
reject = dict(eeg=120e-6, eog=150e-6, grad=200e-12, mag=4e-12)
event_ids = {'cncrt_wrd': 1, 'abs_wrd': 2}
# frequency bands with variable number of cycles for wavelets
stim_delay = 0.034 # delay in s


tmin1=0
tstep1=50
stc_allc=range(17)
vertices_to = [np.arange(10242), np.arange(10242)]
stc_alla=range(17)
stc_alla2=range(17)
stc_allc2=range(17)
Xc_pca=np.zeros((2273724, 17))
Xa_pca=np.zeros((2273724, 17))

vertices_to = [np.arange(10242), np.arange(10242)]
for ii, meg in enumerate(ll):
    print ii
    fname = data_path + meg + 'firstMorphed_SemLoc_icaclean_Concrete_Source_Evoked_m500_700'#'firstMorphed_SemLoc_icaclean_Concrete_Source_Evoked_meanresampled_50_550_100ms_0overlap'
    stcc = mne.read_source_estimate(fname)
    stcc.resample(200)
    stcc.crop(0,0.55)
    fname = data_path + meg + 'firstMorphed_SemLoc_icaclean_Abstract_Source_Evoked_m500_700'#'firstMorphed_SemLoc_icaclean_Abstract_Source_Evoked_meanresampled_50_550_100ms_0overlap' 
    stca = mne.read_source_estimate(fname)
    stca.resample(200)
    stca.crop(0,0.55)
    
    cdata=stcc.data.reshape(stcc.data.shape[0]*stcc.data.shape[1],1)
    adata=stca.data.reshape(stca.data.shape[0]*stca.data.shape[1],1)
    
    Xc_pca[:,ii]=cdata[:,0]
    Xa_pca[:,ii]=adata[:,0]
path_c=data_path+'Xc_pca.mat'  
path_a=data_path+'Xa_pca.mat'    
  
scipy.io.savemat(path_c,{'Xc_pca':Xc_pca})    
scipy.io.savemat(path_a,{'Xa_pca':Xa_pca})    

