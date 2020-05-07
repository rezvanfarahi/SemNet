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
sys.path.insert(1,'/imaging/local/software/mne_python/latest')


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
out_path = '/imaging/rf02/Semnet/stc/permutation/evoked/' # root
uvttest_path = '/imaging/rf02/Semnet/stc/uvttest/evoked/' # root

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

subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14] # removed
#p_inds=[0,1,2,3,4,5,6,7,8,9,10]
p_list=[0.00005,0.05,0.045,0.04,0.03,0.025,0.02,0.01,0.008,0.005,0.002,0.001, 0.0008,0.0005,0.0002,0.0001,0.00005,0.00001]#0.001,0.0009,0.0008,0.0007,0.0006,0.0005]#
print "subject_inds:"
print subject_inds
print "No rejection"

list_all =  ['/meg16_0045/160303/', #8
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
subjects=  [
            'MRI_meg16_0045' ,#8
            'MRI_meg16_0052' ,#10
            'MRI_meg16_0056' ,#11
    	       'MRI_meg16_0069' ,#12
 	       'MRI_meg16_0070' ,#13
            'MRI_meg16_0072' ,#15
            'MRI_meg16_0073' ,#16
            'MRI_meg16_0075' ,#17
            'MRI_meg16_0078' ,#18
            'MRI_meg16_0082' ,#19
            'MRI_meg16_0086' ,#20
            'MRI_meg16_0097' ,#21
            'MRI_meg16_0122' ,#22
            'MRI_meg16_0125' ,#24
            ]

ll = []
for ss in p_inds:
 ll.append(p_list[ss])

print "ll:"
print ll

n_subjects=19
stim_delay = 0.034 # delay in s ((NOTE! not in the example, taken from Olaf's script. Rezvan))
tmin, tmax = -0.1, 0.5
# so we have X: observations * time * vertices in each condition
#X_list=range(2)
event_names = ['colour', 'grey', 'shapes', 'scrambled']#, ,'Pwordc']
#thiscond=0
refcond=3
n_levels=len(event_names)
n_times=151
tmin1=-100
tstep1=4
X=np.zeros((n_subjects,n_times,20484,n_levels))
method='MNE'
for p_threshold in ll:
    
    for ii, meg in enumerate(list_all):
        print ii
        for event_no in range(n_levels):     
            fname_in = data_path + meg + method+'_firstMorphed_ico_newreg_localisers_ica_'+event_names[event_no]+'_Source_Evoked_m100_500'#'_Source_Evoked_m300_600'#'_Source_Evoked_mnrsmpl_50_550_100ms_0ov'
            stc_cond = mne.read_source_estimate(fname_in)
#            stc_cond.resample(40)
#            stc_cond.crop(0.050,0.450)
            X[ii,:,:,event_no]=np.transpose(stc_cond.data,[1,0]) 
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
    for thiscond in np.array([2]):#in range(5):
        Matx1=np.squeeze(X[:,:,:,thiscond]-X[:,:,:,refcond])
        Matx=np.transpose(Matx1,[0,2,1])#np.zeros((n_subjects,20484,n_times))
        MT=np.zeros((20484,n_times))
    
        for cnt in range(Matx.shape[2]):
            print cnt
            MT[:,cnt]= mne.stats.ttest_1samp_no_p(Matx[:,:,cnt], sigma=1e-3, method='relative')
        #MT[np.abs(MT)<t_threshold]=0
    
        vertices_to = [np.arange(10242), np.arange(10242)]
        tval_stc = mne.SourceEstimate(MT, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
        out_file3=uvttest_path + method+'_UVTtest_icomorphed_newreg_14subj_colour_localiser_m100_500_ica_'+event_names[thiscond]+'_'+event_names[refcond]
        tval_stc.save(out_file3)

#col_stc=tval_stc.copy()
shp_stc=tval_stc.copy()
shp_dat=shp_stc.data[:,45:65] 
col_dat=col_stc.data[:,45:65]    
   
st=np.where(np.abs(shp_dat)>3)
ct=np.where(np.abs(col_dat)>3)
tuple_shp=[(st[0][ii],st[1][ii]) for ii in range(len(st[0]))]
tuple_col=[(ct[0][ii],ct[1][ii]) for ii in range(len(ct[0]))]
tuple_ov=[a1 for a1 in tuple_col if a1 in tuple_shp]
list_ov=[c1[0] for c1 in tuple_ov]
final_ov=np.unique(list_ov)

ov_mat=np.zeros((20484,1))
ov_mat[final_ov]=1
ov_stc = mne.SourceEstimate(ov_mat, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
out_file3=uvttest_path + method+'_UVTtest_icomorphed_newreg_14subj_colour_shape_localisers_70_200_overlap'
ov_stc.save(out_file3)



    	