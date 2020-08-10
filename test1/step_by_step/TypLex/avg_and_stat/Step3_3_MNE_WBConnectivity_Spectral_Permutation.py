
"""
=========================================================
TF representation of typlex for frequency bands in source labels
=========================================================

"""
# Authors: Alexandre Gramfort <gramfort@nmr.mgh.harvard.edu>
#
# License: BSD (3-clause)


print(__doc__)

import sys
#sys.path.append('/imaging/local/software/python_packages/nibabel/1.3.0')
sys.path.insert(1,'/imaging/rf02/scikit-learn-0.16.1')
#sys.path.append('/imaging/local/software/python_packages/scikit-learn/v0.14.1')
sys.path.append('/imaging/local/software/mne_python/latest_v0.9')
# for qsub
# (add ! here if needed) 
#sys.path.insert(1,'/imaging/local/software/anaconda/1.9.1/x86_64/envs/py3k')#bin/python')
sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
###
import numpy as np
import mne

import os
import scipy.io

from mne import ( spatial_tris_connectivity, grade_to_tris)
from mne.stats import (spatio_temporal_cluster_1samp_test,
                       summarize_clusters_stc)


###############################################################################
data_path = '/imaging/rf02/TypLexMEG/' # root directory for your MEG data
out_path = '/imaging/rf02/TypLexMEG/icaanalysis_results/typlex/stc/permutation/connectivity/WB_hubs/spectral/newp/'    # where event files are
label_path = '/imaging/rf02/TypLexMEG/fsaverage/label'
sys.path.insert(1,'/imaging/local/software/python_packages/scikit-learn/v0.13.1/lib.linux-x86_64-2.7/sklearn/externals')
import joblib
# 
print sys.argv
p_inds = []
for ss in sys.argv[1:]:
   p_inds.append( int( ss ) )

subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16] # removed
#p_inds=[0,1,2,3,4,5,6,7,8,9,10]
p_list=[0.001,0.0025,0.005,0.008,0.01,0.025,0.05]#[0.001,0.005,0.01,0.025, 0.05]#,0.045,0.04,0.03,0.025,0.02,0.01,0.008,0.005,0.002,0.001,0.0005,0.0001]
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
for ss in p_inds:
 ll.append(p_list[ss])

print "ll:"
print ll



# Labels/ROIs
#labellist = ['atlleft-lh']#, 'atlleftt-rh'
tmin, tmax = -0.5, 0.7
reject = dict(eeg=120e-6, grad=200e-12, mag=4e-12) #eog=150e-6, 
event_ids = {'cncrt_wrd': 1, 'abs_wrd': 2}
# frequency bands with variable number of cycles for wavelets
stim_delay = 0.034 # delay in s
frequencies = np.arange(8, 45, 1)
n_cycles = frequencies / float(7)
n_cycles[frequencies<=20] = 2
n_cycles[frequencies<=20] += np.arange(0,1,0.077)
    # n_cycles[np.logical_and((frequencies>10), (frequencies<=20))] = 3
#'lh.lateralorbitofrontal','rh.lateralorbitofrontal',
labellist=['ATL_latmed-lh','supramarginal-lh','posterior_inferiorfrontal-lh','middle_temporal-lh']
labellistm=['ATLl','SMGl','pIFGl','MTLl']#
#labellist = ['lh.ATL_latmed','lh.supramarginal','lh.posterior_inferiorfrontal','lh.middle_temporal']#[ 'lh.posterior_fusiform','rh.posterior_fusiform', 'lh.lateraloccipital','rh.lateraloccipital','lh.precentral','rh.precentral', 'lh.ATL_latmed','rh.ATL_latmed']#labellist = [ 'leftATL-lh','rightATL-rh']#, 'medialtempright-rh','medialtempleft-lh']'lh.postcentral','rh.postcentral',

#labellist = [ 'lh.posterior_fusiform','rh.posterior_fusiform', 'lh.lateraloccipital','rh.lateraloccipital','lh.precentral','rh.precentral','lh.ATL_latmed','rh.ATL_latmed']#, 'lh.ATLsmall','rh.ATLsmall']#labellist = [ 'leftATL-lh','rightATL-rh']#, 'medialtempright-rh','medialtempleft-lh']'lh.postcentral','rh.postcentral',
#labellist = [ 'lh.posterior_fusiform','rh.posterior_fusiform', 'lh.lateraloccipital','rh.lateraloccipital','lh.precentral','rh.precentral']#labellist = [ 'posterior_fusiform-lh','rh.lateralorbitofron-rh', 'rh.medialorbitofron-rh']#labellist = ['leftATL-rh','rightATL-lh']#, 'rightATL-rh']#'atlleft-lh']#, 'atlright-rh']#,'medialtempright-rh','medialtempleft-lh']
conn_methods=['Coherence','WPLI2d']#,'DWPLI']
bands=['Theta','Alpha','Beta','Gamma']
for p_threshold in ll:
    print p_threshold
    for label_name,lblname in zip(labellist,labellistm):
        print lblname
        for con_method in conn_methods:
            print con_method
            for band in bands:
            	ii=-1
            	#X=np.zeros((20484,2,17))
            	X=np.zeros((20484,2,17))
            	
            	for meg in list_all:
            		ii=ii+1;
            		#print ii
            	
            ### Concrete
            		fname_cncrt = data_path + meg + 'Morphed_ico_typlex_ica_equalized_word_Mlttap_bestV_ctf_'+con_method+'_'+band+'_150_450_200ms_' + label_name[0:-3] 
            		stc_cncrt = mne.read_source_estimate(fname_cncrt)
            		fname_abs = data_path + meg + 'Morphed_ico_typlex_ica_equalized_nonword_Mlttap_bestV_ctf_'+con_method+'_'+band+'_150_450_200ms_' + label_name[0:-3]   
            		stc_abs = mne.read_source_estimate(fname_abs)	
            		stc_data_subtract=np.subtract(np.abs(stc_cncrt.data),np.abs(stc_abs.data))
            		#X[:,:,ii]=stc_subtract.data
            		X[:,:,ii]=stc_data_subtract
            		
            		
            	X = np.transpose(X, [2, 1, 0])	
            	connectivity = spatial_tris_connectivity(grade_to_tris(5))
            	#p_threshold = 0.045
#            	print p_threshold
            	n_subjects=17
            	t_threshold = -scipy.stats.distributions.t.ppf(p_threshold/2., n_subjects - 1)
            	fsave_vertices = [np.arange(10242), np.arange(10242)]
            	n_permutations=1024
            	# I'm replicating the columns because it makes the script faster. Have checked and it doesn't introduce 
            	
#            	print('Computing connectivity.')
            	fname_label = label_path + '/' + 'toremove_wbspokes-lh.label'; labelL = mne.read_label(fname_label)
            	fname_label = label_path + '/' + 'toremove_wbspokes-rh.label'; labelR = mne.read_label(fname_label)
            	labelss=labelL+labelR
            	bb=stc_cncrt.in_label(labelss)
            	fsave_vertices = [np.arange(10242), np.arange(10242)]
            	nnl=np.in1d(fsave_vertices[0],bb.lh_vertno)
            	nnr=np.in1d(fsave_vertices[1],bb.rh_vertno)
            	spatial_exclude=np.hstack((fsave_vertices[0][nnl], fsave_vertices[0][nnr]+10242))
                
            
            	
            	print('Clustering ' +band); max_step=1;
            	T_obst, clusterst, cluster_p_valuest, H0t = clut = spatio_temporal_cluster_1samp_test(X, connectivity=connectivity,  n_permutations=n_permutations, n_jobs=4, threshold=t_threshold, tail=0,t_power=1, spatial_exclude=spatial_exclude,step_down_p=0.05)#,max_step=0)
            	#    Now select the clusters that are sig. at p < 0.05 (note that this value
            	#    is multiple-comparisons corrected).
            	clus_p_values=(p_threshold/0.05)*cluster_p_valuest
               
            	good_cluster_indst = np.where(clus_p_values <= 0.003125)[0]
            	print cluster_p_valuest.min(), clus_p_values.min()
            	
            	
            	
            
            	tstep=1e-3
            	fsave_vertices = [np.arange(10242), np.arange(10242)]
            	if good_cluster_indst.size>0:
                	p_thresh=np.max(cluster_p_valuest[good_cluster_indst])#cluster_p_valuest.min()+0.000001
            	
            	if clus_p_values.min()<=0.003125:
                    print p_thresh
                    
                    stc_all_cluster_vist = summarize_clusters_stc(clut, tstep=tstep, p_thresh=p_thresh+0.000001, vertno=fsave_vertices, subject='fsaverage')
                    out_file1=out_path + 'ClusP_Morphed_ico_TL_ica_eq_sub_Mlttap_bestV_ctf_'+con_method+'_'+band+'_150_450_200ms_sx_ms' + str(max_step) +'_'+ str(p_threshold)[2:5] + '_' + str(p_thresh)[2:5]+ '_' + lblname
                    stc_all_cluster_vist.save(out_file1)
                      	
                    Matx=np.zeros((20484,2))
                    T=np.divide(T_obst,np.absolute(T_obst))
                    for cc in range(good_cluster_indst.shape[0]):
                        Matx[clusterst[good_cluster_indst[cc]][1],clusterst[good_cluster_indst[cc]][0]]=T[clusterst[good_cluster_indst[cc]][0],clusterst[good_cluster_indst[cc]][1]]
                    
                    tmin1=150
                    tstep1=100
                    vertices_to = [np.arange(10242), np.arange(10242)]
                    matx_stc = mne.SourceEstimate(Matx, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
                    out_file2=out_path + 'ClusP_sw_Morphed_ico_TL_ica_eq_sub_Mlttap_bestV_ctf_'+con_method+'_'+band+'_150_450_200ms_sx_ms' + str(max_step) +'_'+ str(p_threshold)[2:5] + '_' + str(p_thresh)[2:5]+ '_' + lblname
                    matx_stc.save(out_file2)
            	
