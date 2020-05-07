"""
=================================================================
Permutation t-test on source data with spatio-temporal clustering
=================================================================

Tests if the evoked response is significantly different between
conditions across subjects (simulated here using one subject's data).
The multiple comparisons problem is addressed with a cluster-level
permutation test across space and time.

"""

# Authors: Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
#          Eric Larson <larson.eric.d@gmail.com>
# License: BSD (3-clause)

print(__doc__)
# Russell's addition
import sys
#sys.path.insert(1,'/imaging/rf02/Semnet/mne_python0.9')
#sys.path.insert(1,'/imaging/local/software/mne_python/latest')
sys.path.insert(1,'/imaging/local/software/mne_python/latest_v0.11')

import os.path as op
import numpy as np
from numpy.random import randn
from scipy import stats as stats

import mne
from mne import (io, spatial_tris_connectivity, compute_morph_matrix,spatio_temporal_tris_connectivity,
                 grade_to_tris)
from mne.epochs import equalize_epoch_counts
from mne.stats import (spatio_temporal_cluster_1samp_test,
                       summarize_clusters_stc)
from mne.minimum_norm import apply_inverse, read_inverse_operator
from mne.stats import f_threshold_mway_rm, f_mway_rm, fdr_correction, spatio_temporal_cluster_test
from mne.datasets import sample
from mne.viz import mne_analyze_colormap
import os

###############################################################################
# Set parameters
out_path = '/imaging/rf02/Semnet/semnet4semloc/stc/permutation/evoked/' # root
uvttest_path = '/imaging/rf02/Semnet/semnet4semloc/stc/uvttest/evoked/1_48/' # root

if not os.path.exists(out_path):
    os.makedirs(out_path)
if not os.path.exists(uvttest_path):
    os.makedirs(uvttest_path)
data_path = '/imaging/rf02/Semnet/'	# where subdirs for MEG data are
inv_path = '/imaging/rf02/Semnet/'
print sys.argv
p_inds = [0]
for ss in sys.argv[1:]:
   p_inds.append( int( ss ) )
label_path = '/imaging/rf02/TypLexMEG/fsaverage/label'

subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19] # removed
#p_inds=[0,1,2,3,4,5,6,7,8,9,10]
p_list=[0.00005,0.05,0.045,0.04,0.03,0.025,0.02,0.01,0.008,0.005,0.002,0.001, 0.0008,0.0005,0.0002,0.0001,0.00005,0.00001]#0.001,0.0009,0.0008,0.0007,0.0006,0.0005]#
print "subject_inds:"
print subject_inds
print "No rejection"
b='gamma'
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
for ss in p_inds:
 ll.append(p_list[ss])

print "ll:"
print ll

n_subjects=19
stim_delay = 0.034 # delay in s ((NOTE! not in the example, taken from Olaf's script. Rezvan))
tmin, tmax = -0.3, 0.6
# so we have X: observations * time * vertices in each condition
#X_list=range(2)
event_names = ['Neutral', 'Emotional','Concrete']#, ,'Pwordc']'Neutral', 'Emotional',
#thiscond=0
thiscond=2
n_levels=len(event_names)
n_times=901
tmin1=50
tstep1=100
X=np.zeros((n_subjects,n_times,20484,n_levels))
nwins=5
Xmean=np.zeros((n_subjects,nwins,20484,n_levels))
for p_threshold in ll:
    
    for ii, meg in enumerate(list_all):
        print ii
        for event_no in range(n_levels):     
            fname_in = data_path + meg + 'firstMorphed_ico_oldreg_SemDec_SL_1_48ica_'+event_names[event_no]+'_Source_Evoked_m300_600'#_'+b#'firstMorphed_ico_SemDec_ica_'+event_names[event_no]+'_Source_Evoked_mnrsmpl_50_550_100ms_0ov'#'_Source_Evoked_m300_600'#'_Source_Evoked_mnrsmpl_50_550_100ms_0ov'
            stc_cond = mne.read_source_estimate(fname_in)
#            stc_cond.resample(100)
#            stc_cond.crop(0.050,0.450)
            X[ii,:,:,event_no]=np.transpose(stc_cond.data,[1,0]) #[:,350:650]
    X_list=[np.squeeze(x) for x in np.split(X, n_levels, axis=-1)]
    factor_levels = [n_levels]  # number of levels in each factor
    effects = 'A'  # this is the default signature for computing all effects
    return_pvals = False          
    def stat_fun(*args):  
        return f_mway_rm(np.swapaxes(args, 1, 0), factor_levels=factor_levels,effects=effects, return_pvals=return_pvals)[0]

#            xx=np.asarray(args)
#            xx1=np.reshape(xx,(xx.shape[0],xx.shape[1],xx.shape[2]*xx.shape[3]))
#            xx2=np.transpose(xx1,[1,0,2])
#            return f_mway_rm(xx2, factor_levels=factor_levels,effects=effects, return_pvals=False)[0]
    print('Visualizing clusters.')
    wcnt=-1
    for wcnt1,wcnt2 in zip([350,450,550,650,750],[450,550,650,750,850]):#range(nwins):
        wcnt=wcnt+1
        Xmean[:,wcnt,:,:]=np.mean(X[:,wcnt1:wcnt2,:,:],1)
    for refcond in np.array([0,1]):#in range(5):
        Matx1=np.squeeze(np.subtract(Xmean[:,:,:,thiscond],Xmean[:,:,:,refcond]))
        Matx=np.transpose(Matx1,[0,2,1])#np.zeros((n_subjects,20484,n_times))
        MT=np.zeros((20484,nwins))
    
        for cnt in range(Matx.shape[2]):
            print cnt
            MT[:,cnt]= mne.stats.ttest_1samp_no_p(Matx[:,:,cnt], sigma=1e-3, method='relative')
        #MT[np.abs(MT)<t_threshold]=0
    
        vertices_to = [np.arange(10242), np.arange(10242)]
        tval_stc = mne.SourceEstimate(MT, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
        out_file3=uvttest_path + 'UVTtest_icomorphed_oldreg_19subj_SemDec_m300_600_100ms_SL_1_48ica_'+event_names[thiscond]+'_'+event_names[refcond]#+'_'+b
        tval_stc.save(out_file3)
    
    	