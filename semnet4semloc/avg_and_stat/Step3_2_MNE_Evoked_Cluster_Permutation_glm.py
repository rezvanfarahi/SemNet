"""
=================================================================
Permutation rm ANOVA on source data with spatio-temporal clustering
=================================================================

Tests if the evoked response is significantly different between
conditions across subjects (simulated here using one subject's data).
The multiple comparisons problem is addressed with a cluster-level
permutation test across space and time.

"""

# Authors: Rezvan Farahibozorg, July 2020

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
from nipy.modalities.fmri.glm import GeneralLinearModel


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
import copy

###############################################################################
# Set parameters
exclude_wbmedial=True
exclude_ROIs=False
win_20ms=False
if exclude_wbmedial:
    out_path = '/imaging/rf02/Semnet/semnet4semloc/stc/permutation/evoked/glm/' # root
if exclude_ROIs:
    out_path = '/imaging/rf02/Semnet/semnet4semloc/stc/permutation/masked_ROIs/evoked/glm/' # root
uvttest_path = '/imaging/rf02/Semnet/semnet4semloc/stc/uvttest/evoked/glm/' # root
if win_20ms:
    out_path=out_path+'win_20ms/'
    uvttest_path=uvttest_path+'win_20ms/'
if not os.path.exists(out_path):
    os.makedirs(out_path)
if not os.path.exists(uvttest_path):
    os.makedirs(uvttest_path)
data_path = '/imaging/rf02/'#Semnet/'	# where subdirs for MEG data are
inv_path = '/imaging/rf02/'#Semnet/'
print (sys.argv)
p_inds = [0]
for ss in sys.argv[1:]:
   p_inds.append( int( ss ) )
label_path = '/imaging/rf02/TypLexMEG/fsaverage/label'

#subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19] # removed
#p_inds=[0,1,2,3,4,5,6,7,8,9,10]
p_list=[0.005]#,0.045,0.04,0.03,0.025,0.02,0.01,0.008,0.005,0.002,0.001, 0.0008,0.0005,0.0002,0.0001,0.00005,0.00001]#0.001,0.0009,0.0008,0.0007,0.0006,0.0005]#
#print ("subject_inds:")
#print (subject_inds)
print ("No rejection")
list_all={'semloc' :  ['/TypLexMEG/meg10_0378/101209/', 
            '/TypLexMEG/meg10_0390/101214/',
            '/TypLexMEG/meg11_0026/110223/', 
            '/TypLexMEG/meg11_0050/110307/', 
            '/TypLexMEG/meg11_0052/110307/', 
            '/TypLexMEG/meg11_0069/110315/', 
            '/TypLexMEG/meg11_0086/110322/', 
            '/TypLexMEG/meg11_0091/110328/', 
            '/TypLexMEG/meg11_0096/110404/', 
            '/TypLexMEG/meg11_0101/110411/', 
            '/TypLexMEG/meg11_0102/110411/', 
            '/TypLexMEG/meg11_0112/110505/',
            '/TypLexMEG/meg11_0104/110412/',  
            '/TypLexMEG/meg11_0118/110509/', 
            '/TypLexMEG/meg11_0131/110519/', 
            '/TypLexMEG/meg11_0144/110602/', 
            '/TypLexMEG/meg11_0147/110603/', 
            ],
'semnet1' :  ['/Semnet/meg16_0032/160218/', #1  '/meg16_0030/160216/', #0
            '/Semnet/meg16_0034/160219/', #3
            '/Semnet/meg16_0035/160222/', #4
            '/Semnet/meg16_0042/160229/', #7
            '/Semnet/meg16_0045/160303/', #8
            '/Semnet/meg16_0052/160310/', #10
            '/Semnet/meg16_0056/160314/',#11
            '/Semnet/meg16_0069/160405/',#12 
            '/Semnet/meg16_0070/160407/', #13
            '/Semnet/meg16_0072/160408/', #15
            '/Semnet/meg16_0073/160411/', #16
            '/Semnet/meg16_0075/160411/', #17
            '/Semnet/meg16_0078/160414/', #18
            '/Semnet/meg16_0082/160418/', #19
            '/Semnet/meg16_0086/160422/', #20
            '/Semnet/meg16_0097/160512/', #21 
            '/Semnet/meg16_0122/160707/', #22 
            '/Semnet/meg16_0125/160712/', #24 
            ],

'semnet2' :  ['/Semnet/meg16_0032/160218/', #1  '/meg16_0030/160216/', #0
            '/Semnet/meg16_0034/160219/', #3
            '/Semnet/meg16_0035/160222/', #4
            '/Semnet/meg16_0042/160229/', #7
            '/Semnet/meg16_0045/160303/', #8
            '/Semnet/meg16_0052/160310/', #10
            '/Semnet/meg16_0056/160314/',#11
            '/Semnet/meg16_0069/160405/',#12 
            '/Semnet/meg16_0070/160407/', #13
            '/Semnet/meg16_0072/160408/', #15
            '/Semnet/meg16_0073/160411/', #16
            '/Semnet/meg16_0075/160411/', #17
            '/Semnet/meg16_0078/160414/', #18
            '/Semnet/meg16_0082/160418/', #19
            '/Semnet/meg16_0086/160422/', #20
            '/Semnet/meg16_0097/160512/', #21 
            '/Semnet/meg16_0122/160707/', #22 
            '/Semnet/meg16_0125/160712/', #24 
            ]
}



ll = []
for ss in p_inds:
 ll.append(p_list[ss])

print ("ll:")
print (ll)

#n_subjects=len(subjects)
# so we have X: observations * time * vertices in each condition
#X_list=range(2)
semtasks=['SemDec','LD']#
event_names = {'semloc':['Concrete','Abstract'],'semnet1':['Concrete','Emotional'],'semnet2':['Concrete','Emotional']}#, ,'Pwordc']'Neutral', 'Emotional',
effect_names=['contrast','interaction']#,'task']
#all_effects=['B','A:B','A']

vertices_avg = [np.arange(10242), np.arange(10242)]
n_levels=len(semtasks)
all_nconds=[len(event_names['semloc']),len(event_names['semnet1']),len(event_names['semnet2'])]
nconds=copy.deepcopy(all_nconds[0])
factor_levels = [n_levels,n_levels]  # number of levels in each factor
grand_tstep=100
grand_tmin=350
grand_tmax=751-grand_tstep
n_times=len(list(np.arange(350,651,grand_tstep)))
tmin1=50
tstep1=copy.deepcopy(grand_tstep)
Nv = len(vertices_avg[0])+len(vertices_avg[1])
Nt=copy.deepcopy(n_times)
Ns = [len(list_all['semloc']),len(list_all['semnet1']),len(list_all['semnet2'])]
Nsubj=Ns[0]+Ns[1]+Ns[2]

X=np.zeros((Nsubj,nconds,Nv,Nt))#np.zeros((n_subjects,n_times,20484,n_levels))
#Xmean=np.zeros((n_subjects,nwins,20484,n_levels))
for p_threshold in ll: 
    ii=-1  
    for taski,task_name in enumerate(['semloc','semnet1','semnet2']):
        for meg in list_all[task_name]:
            ii=ii+1
            tecnt=-1
            for evcnt, event_name in enumerate(event_names[task_name]):#contrast is B
                tecnt=tecnt+1
                if task_name=='semloc':
                    fname_in = data_path + meg + 'firstMorphed_ico_SemLoc_ica_'+event_name+'_Source_Evoked_m500_700'
                if task_name=='semnet1':
                    fname_in = data_path + meg + 'firstMorphed_ico_oldreg_SemDec_SL_1_48ica_'+event_name+'_Source_Evoked_m300_600'
                if task_name=='semnet2':
                    fname_in = data_path + meg + 'firstMorphed_ico_oldreg_LD_SL_1_48ica_'+event_name+'_Source_Evoked_m300_600'
                
                stc_cond = mne.read_source_estimate(fname_in)
                #            stc_cond.resample(100)
                #            stc_cond.crop(0.050,0.450)
                if task_name=='semloc':
                    wcnt=-1
                    for wcnt1,wcnt2 in zip(list(np.arange(grand_tmin+200,grand_tmax+200,grand_tstep)),list(np.arange(grand_tmin+grand_tstep+200,grand_tmax+grand_tstep+200,grand_tstep))):#range(nwins):
                        print (wcnt1,wcnt2)
                        wcnt=wcnt+1
                        X[ii,evcnt,:,wcnt]=np.mean(stc_cond.data[:,wcnt1:wcnt2],1)
                    else:
                        wcnt=-1
                        for wcnt1,wcnt2 in zip(list(np.arange(grand_tmin,grand_tmax,grand_tstep)),list(np.arange(grand_tmin+grand_tstep,grand_tmax+grand_tstep,grand_tstep))):#range(nwins):
                            print (wcnt1,wcnt2)
                            wcnt=wcnt+1
                            X[ii,evcnt,:,wcnt]=np.mean(stc_cond.data[:,wcnt1:wcnt2],1)

                #X[ii,:,:,event_no]=np.transpose(stc_cond.data,[1,0]) #[:,350:650]
    X1=np.transpose(X, [0, 3, 2, 1]).copy()#X1 needs to be (Nsub, Ntime, Nvox, Ncond)
    X_list=[np.squeeze(x) for x in np.split(X1, nconds, axis=-1)]#Xlist is now (Ncond,Nsub, Ntime, Nvox)
    

    source_space = grade_to_tris(5)
    print('Computing connectivity.')
    connectivity = spatial_tris_connectivity(source_space)

    pthresh = p_threshold
    # max_step=1;
    

    #    To speed things up a bit we will ...
    n_permutations = 5000  # ... run fewer permutations (reduces sensitivity)
    fsave_vertices = [np.arange(10242), np.arange(10242)]
    if exclude_wbmedial:
        fname_label = label_path + '/' + 'toremove_wbspokes-lh.label'; labelL = mne.read_label(fname_label)
        fname_label = label_path + '/' + 'toremove_wbspokes-rh.label'; labelR = mne.read_label(fname_label)       
    if exclude_ROIs:
        fname_label='/imaging/rf02/Semnet/semnet4semloc//mask_labels_ATL_IFG_MTG_AG-lh.label'; labelL = mne.read_label(fname_label)
        fname_label='/imaging/rf02/Semnet/semnet4semloc//mask_labels_ATL_IFG_MTG_AG-rh.label'; labelR = mne.read_label(fname_label)
        #labelmask=mne.read_label(fname_label,subject='fsaverage')
        #labelmask.values.fill(1.0)
        #bb=stc_cond.in_label(labelmask)
        #nnl=np.in1d(fsave_vertices[0],bb.lh_vertno)
        #spatial_exclude=fsave_vertices[0][nnl].copy()
    labelss=labelL+labelR
    bb=stc_cond.in_label(labelss)  
    nnl=np.in1d(fsave_vertices[0],bb.lh_vertno)
    nnr=np.in1d(fsave_vertices[1],bb.rh_vertno)
    spatial_exclude=np.hstack((fsave_vertices[0][nnl], fsave_vertices[0][nnr]+10242))
    print('Clustering.')
    #    t_threshold = -stats.distributions.t.ppf(p_threshold/2., n_subjects - 1)#dict(start=0, step=.1)#
    #t_threshold=2
    tail=0
    max_step=5
    desing_mat = np.hstack((np.ones((Ns[0],1)),np.zeros((Ns[0],2))))
    desing_mat = np.vstack((desing_mat, np.hstack((np.zeros((2*Ns[1],1)), np.kron(np.eye(2),np.ones((Ns[1],1)))))))
    desing_mat = np.hstack((desing_mat ,np.vstack((np.zeros((Ns[0],Ns[1])), np.eye(Ns[1]), np.eye(Ns[1])))))

    effects={'contrast':np.hstack((np.asarray([1, 1, 1])[np.newaxis,:],np.ones((1,Ns[1]))*2/Ns[1] )),#this is for the main effect
    'interaction':np.hstack((np.eye(3)-np.mean(np.eye(3)),np.zeros((3,Ns[1]))))}    
    for effect_name in effect_names:
        f_thresh = 2 #f_threshold_mway_rm(n_subjects, factor_levels, this_effect, pthresh)
        #NOTEE! stat_fun will have to return a 1-D array 
        def stat_fun(*args): #this function swaps Ncond and Nsub in X_list; that is the input dimension that anova requires
            cond1=args[0].copy()#cond1 is Nsub x Ntime x Nv
            cond1=cond1.reshape(cond1.shape[0], np.prod(cond1.shape[1:]))#now cond1 is Nsub x Ntime*Nv
            cond2=args[1].copy()
            cond2=cond2.reshape(cond2.shape[0], np.prod(cond2.shape[1:]))#now cond2 is Nsub x Ntime*Nv
            cond_subtract=cond1-cond2
            model = GeneralLinearModel(desing_mat)
            model.fit(cond_subtract)
            glm_pvals = model.contrast(effects[effect_name]).p_value()
            glm_tvalues = model.contrast(effects[effect_name]).stat()
            return glm_tvalues

        #effects = 'B'  # A*B is the default signature for computing all effects, A here is task effect, B contrast 
        return_pvals = False         
        T_obs, clusters, cluster_p_values, H0 = clu = \
        spatio_temporal_cluster_test(X_list, connectivity=connectivity, n_jobs=4,step_down_p=0.05,max_step=max_step,t_power=1,
                                        threshold=f_thresh, stat_fun=stat_fun, spatial_exclude=spatial_exclude,
                                        n_permutations=n_permutations,tail=0,
                                        buffer_size=None)
    
        print('Visualizing clusters.')
        if len(cluster_p_values)>-1:
            if cluster_p_values.min()<1:
                p_thr=cluster_p_values.min()
                good_cluster_inds = np.where(cluster_p_values <=p_thr)[0]
                print (cluster_p_values[good_cluster_inds]); print (good_cluster_inds)
                stc_all_cluster_vis = summarize_clusters_stc(clu, tstep=1e-3 * tstep1, vertices=fsave_vertices, subject='fsaverage', p_thresh=p_thr+0.0001)
                
                out_file1=out_path + 'ClusPer_GLM_Evoked_icomorphed_oldreg_clusterp'+str(p_threshold)[2:]+'_p'+str(p_thr)[2:]+'_53subj_pnt1_48ica_'+effect_name+'_'+str(max_step)
                stc_all_cluster_vis.save(out_file1)
                
                Matx=np.zeros((20484,n_times))
                T=np.divide(T_obs,np.absolute(T_obs))
                for cc in range(good_cluster_inds.shape[0]):
                    Matx[clusters[good_cluster_inds[cc]][1],clusters[good_cluster_inds[cc]][0]]=T[clusters[good_cluster_inds[cc]][0],clusters[good_cluster_inds[cc]][1]]
                    
                
                #aa=Matx[clusters[good_cluster_inds[cc]][1],clusters[good_cluster_inds[cc]][0]]
                #bb=T_obs[clusters[good_cluster_inds[cc]][0],clusters[good_cluster_inds[cc]][1]]<0
                #aa[bb]=-1
                #Matx[clusters[good_cluster_inds[cc]][1],clusters[good_cluster_inds[cc]][0]]=aa
                
                
                matx_stc = mne.SourceEstimate(Matx, vertices=vertices_avg,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
                out_file2=out_path + 'ClusPer_GLM_Evoked_sw_icomorphed_oldreg_clusterp'+str(p_threshold)[2:]+'_p'+str(p_thr)[2:]+'_53subj_pnt1_48ica_'+effect_name+'_'+str(max_step)
                matx_stc.save(out_file2)
        

    	
