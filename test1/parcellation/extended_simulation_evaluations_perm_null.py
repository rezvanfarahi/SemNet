# -*- coding: utf-8 -*-
"""
Created on Mon May 15 01:19:26 2017

@author: rf02
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Dec 23 14:34:39 2015

@author: rf02
"""

"""
==========================================================
Compute point-spread functions (PSFs) for MNE/dSPM/sLORETA
==========================================================

PSFs are computed for four labels in the MNE sample data set
for linear inverse operators (MNE, dSPM, sLORETA).
PSFs describe the spread of activation from one label
across the cortical surface.
"""

# Authors: Olaf Hauk <olaf.hauk@mrc-cbu.cam.ac.uk>
#          Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
#
# License: BSD (3-clause)



# Author: Martin Luessi <mluessi@nmr.mgh.harvard.edu>
#
# License: BSD (3-clause)

print(__doc__)
print "hi! this script is for spectral functional connectivity. Very nice maps are fruits which worth waiting for I bet!"
import sys

#sys.path.append('/imaging/local/software/python_packages/nibabel/1.3.0')
#sys.path.append('/imaging/local/software/python_packages/pysurfer/v0.4')
#sys.path.append('/imaging/rf02/scikit-learn-0.15.0')
#
##sys.path.append('/imaging/local/software/python_packages/scikit-learn/v0.14.1')
#sys.path.append('/imaging/local/software/mne_python/la_v0.9')
#sys.path.insert(1,'/imaging/rf02/mne_python_11')
#sys.path.insert(1,'/imaging/local/software/mne_python/la_v0.9full')
sys.path.insert(1,'/imaging/local/software/mne_python/latest_v0.13')


sys.path.insert(1,'/imaging/rf02/scikit-learn-0.16.1')
#sys.path.insert(1,'/home/rf02/.local/lib/python2.7/site-packages')


# for qsub
# (add ! here if needed) 
#sys.path.insert(1,'/imaging/local/software/EPD/latest/x86_64/bin/python')
#sys.path.insert(1,'/imaging/local/software/anaconda/2.4.1/2/lib/python2.7/site-packages')
#sys.path.insert(1,'/imaging/local/software/anaconda/latest/x86_64/pkgs/python-2.7.13-0')
#sys.path.insert(1,'/imaging/local/software/anaconda/latest/x86_64')
#sys.path.insert(1,'/imaging/local/software/anaconda/latest/3/lib/python3.5/site-packages')
sys.path.insert(1,'/imaging/local/software/anaconda/2.4.1/2/lib/python2.7')
sys.path.insert(1,'/imaging/local/software/python_packages/scikit-learn/v0.13.1/lib.linux-x86_64-2.7/sklearn/externals')
import joblib
#sys.path.append('/imaging/rf02/TypLexMEG/mne_python_v8/')
#sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
###
#import sklearn
#from scipy.stats.mstats import zscore
import os
import numpy as np
#reload(np)
print np.__path__
print np.__version__
import mne
reload(mne)
#print mne.__path__
from mne import io
#from mne.io import Raw
#import operator
from mne.minimum_norm import (read_inverse_operator, apply_inverse_epochs)
import scipy.io
from mne.label import _n_colors as nc
from mne.label import _read_annot as mnera
from mne.label import _write_annot as mnewa
from mne.connectivity import spectral_connectivity
from os import walk

import matplotlib.pyplot as plt
import copy
import pickle
print mne.__path__
import time
import math
from scipy import stats as scistats

#def make_parcellation_mine(src,labels,colors,subject,fname):
#    annotl=np.zeros((src[0]['np'],1))
#    annotr=np.zeros((src[1]['np'],1))
#    CTAB_L=np.zeros((len(labels)/2,5))
#    CTAB_R=np.zeros((len(labels)/2,5))
#    annot_id_coding = np.array((1, 2 ** 8, 2 ** 16))
#    annot_id=np.sum(colors[:,:3] * annot_id_coding, axis=1)
#    for lcntn in range(len(labels)):
#        print lcntn
#        if labels[lcntn].name.endswith('lh'):
#            annotl[labels[lcntn].vertices]=annot_id[lcntn/2]
#        else:
#            annotr[labels[lcntn].vertices]=annot_id[int(lcntn/2)]
#    CTAB_L[:,:4]=colors; CTAB_L[:,4]=annot_id
#    CTAB_R[:,:4]=colors; CTAB_R[:,4]=annot_id
#    label_names=[labels[lh].name[:-3] for lh in range(len(labels))]
#    NAMES_L=label_names[::2]
#    NAMES_R=label_names[1::2]
#    ## write morphed labels
#    fnamel=data_path + subject+'/label/lh.'+subject+'_'+fname+'.annot'
#    mnewa(fnamel, np.squeeze(annotl), CTAB_L, NAMES_L)
#    
#    fnamer=data_path + subject+'/label/rh.'+subject+'_'+fname+'.annot'
#    mnewa(fnamer, np.squeeze(annotr), CTAB_R, NAMES_R)

data_path = '/imaging/rf02/TypLexMEG/' # root directory for your MEG data
os.chdir(data_path)
event_path = '/group/erp/data/olaf.hauk/MEG/TypLexM/data/MF22_new/'
#label_path = '/group/erp/data/olaf.hauk/MEG/TypLexM/labels/createdlabels_SL/'
inv_path = '/imaging/rf02/TypLexMEG/'
subjects_dir = data_path 
print sys.argv
perm_inds = []
for ss in sys.argv[1:]:
   perm_inds.append( int( ss ) )

subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
print "subject_inds:"
print subject_inds
print "No rejection"
t0=time.time()
n_sims=70
sim_offset=0
sim_all=list(np.asarray(range(n_sims))+sim_offset)
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
'meg11_0147/110603/'
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
for ss in perm_inds:
 ll.append(sim_all[ss])


print "ll:"
print ll

label_path='/imaging/rf02/TypLexMEG/parcellation/labels_and_CTFs/'
results_path='/imaging/rf02/parcellation/leakage_pattern/'

##grow ROIs
srcin='/imaging/rf02/TypLexMEG/fsaverage/bem/fsaverage_dist-ico-5-src.fif'
#src_avg2=src_avg.copy()
src_avg = mne.read_source_spaces(srcin)
subject_avg='fsaverage'
#nseeds=20
ext=25
#connperc=0.25#0.25


#labellist_path = [label_path+'rostfront_smallseed-lh.label',label_path+'midtemp_smallseed-lh.label',label_path+'latocci_smallseed-lh.label']#[ 'lh.posterior_fusiform','rh.posterior_fusiform', 'lh.lateraloccipital','rh.lateraloccipital','lh.precentral','rh.precentral', 'lh.ATL_latmed','rh.ATL_latmed']#labellist = [ 'leftATL-lh','rightATL-rh']#, 'medialtempright-rh','medialtempleft-lh']'lh.postcentral','rh.postcentral',
#labellist_path = [label_path+'IFG_smallseed-lh.label',label_path+'IPL_smallseed-lh.label']#[ 'lh.posterior_fusiform','rh.posterior_fusiform', 'lh.lateraloccipital','rh.lateraloccipital','lh.precentral','rh.precentral', 'lh.ATL_latmed','rh.ATL_latmed']#labellist = [ 'leftATL-lh','rightATL-rh']#, 'medialtempright-rh','medialtempleft-lh']'lh.postcentral','rh.postcentral',
#labellist_path = [label_path+'LO_smallseed-lh.label',label_path+'RMF_smallseed-lh.label',label_path+'STS_smallseed-lh.label']#[ 'lh.posterior_fusiform','rh.posterior_fusiform', 'lh.lateraloccipital','rh.lateraloccipital','lh.precentral','rh.precentral', 'lh.ATL_latmed','rh.ATL_latmed']#labellist = [ 'leftATL-lh','rightATL-rh']#, 'medialtempright-rh','medialtempleft-lh']'lh.postcentral','rh.postcentral',
#labellist=[mne.read_label(label) for label in labellist_path]
#for lbl in labellist_m:
#    lbl.values.fill(1.0)
##
#which_parc=1 #names below
    

            
            
            
t=time.time()
            
#labels_pre.remove(labels_pre[-1])
#    labels_pre.remove(labels_pre[-1])
TP_grand=np.zeros((4,3,5,9));#the 9: COH:leak3,noise3,leak1,noise1, imCOH:leak3,noise3,ideal,leak1,noise1,
FP_grand=np.zeros((4,3,5,9));
TN_grand=np.zeros((4,3,5,9));
FN_grand=np.zeros((4,3,5,9));
seedcnt=-1;
for nseeds in [3,5,10,15]:
    print nseeds
    seedcnt=seedcnt+1;
    concnt=-1;
    for nconns in [1,0.5,0.25]:
        print nconns
        concnt=concnt+1;
        for which_parc in range(5):
            nsims=35;#length(FileName);
            nsubs=17;
            print which_parc
            #%both SNRs
            if which_parc==0:
                    nnodes=68
                    inpath='/imaging/rf02/parcellation/extended_simulations/sims70_offset0_seeds'+str(nseeds)+'_conns'+str(nconns)+'_bothsnrs_ext25/aparc/'
                    filelist0=list()
                    for(dirpath,dirnames,filenames) in walk(inpath):
                        filelist0.extend(filenames)
                    filelist=[fn for fn in filelist0 if fn[:12]=='parc_connmat']
                            
#        for filename in filenames:
#            if not filename.endswith('.mat'):
#                fcnt=fcnt+1
#                filepath=in_path+filename
#                with open (filepath, 'rb') as fp:
#                    itemlist = pickle.load(fp)
#                thisvect=np.asarray(itemlist)
#                print fcnt, filename
#                missedcons[fcnt,:]=thisvect.copy()
#                    cd(inpath)
#                    aa=dir ('*parc_connmat*.mat');
            elif which_parc==1:                    
                    nnodes=74;
                    inpath='/imaging/rf02/parcellation/extended_simulations/sims70_offset0_seeds'+str(nseeds)+'_conns'+str(nconns)+'_bothsnrs_ext25/aparc_mod/'
                    filelist0=list()
                    for(dirpath,dirnames,filenames) in walk(inpath):
                        filelist0.extend(filenames)
                    filelist=[fn for fn in filelist0 if fn[:12]=='parc_connmat']                   
            elif which_parc==2:
                    nnodes=148;
                    inpath='/imaging/rf02/parcellation/extended_simulations/sims70_offset0_seeds'+str(nseeds)+'_conns'+str(nconns)+'_bothsnrs_ext25/aparc2009/'
                    filelist0=list()
                    for(dirpath,dirnames,filenames) in walk(inpath):
                        filelist0.extend(filenames)
                    filelist=[fn for fn in filelist0 if fn[:12]=='parc_connmat']
            elif which_parc==3:
                    nnodes=74;
                    inpath='/imaging/rf02/parcellation/extended_simulations/sims70_offset0_seeds'+str(nseeds)+'_conns'+str(nconns)+'_bothsnrs_ext25/aparc2009_mod/'
                    filelist0=list()
                    for(dirpath,dirnames,filenames) in walk(inpath):
                        filelist0.extend(filenames)
                    filelist=[fn for fn in filelist0 if fn[:12]=='parc_connmat']
            elif which_parc==4:
                    nnodes=70;
                    inpath='/imaging/rf02/parcellation/extended_simulations/sims70_offset0_seeds'+str(nseeds)+'_conns'+str(nconns)+'_bothsnrs_ext25/RG/'
                    filelist0=list()
                    for(dirpath,dirnames,filenames) in walk(inpath):
                        filelist0.extend(filenames)
                    filelist=[fn for fn in filelist0 if fn[:12]=='parc_connmat']
                        
            missed_path='/imaging/rf02/parcellation/extended_simulations/sims70_offset0_seeds'+str(nseeds)+'_conns'+str(nconns)+'_bothsnrs_ext25/missedcons/allmissedcons.mat'
            missedcons_pre=scipy.io.loadmat(missed_path,mat_dtype=True); missedcons=np.squeeze(missedcons_pre['missedcons'])
#            this_pathto=out_path+outname[0]+'parc_connmat_sim'+str(perm_inds[0]+sim_offset)+'.mat'#'DKA68_orig_modified_coh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot            
#        scipy.io.savemat(this_pathto,{'parc0_connmat':parc0_connmat})
#            missedcons=np.squeeze(np.mean(missedcons,axis=2));
            missedcons=missedcons[:,which_parc]
            connmat=np.zeros((nnodes,nnodes,10,nsubs,nsims))#;%leakageSNR3,noiseSNR3,ideal,leakageSNR1,noiseSNR1 (COH 1st, imCOH2nd)
#            cc=-1;
            for fcnt,filename in enumerate(filelist):#=1:nsims:
#                cc=cc+1
                filepath=inpath+filename
                conpre=scipy.io.loadmat(filepath,mat_dtype=True); 
                connmat[:,:,:,:,fcnt]=np.squeeze(conpre['parc'+str(which_parc)+'_connmat'])

            if which_parc==0:
                inpath2='/imaging/rf02/parcellation/extended_simulations/null_sims40_offset0_bothsnrs/aparc/'
                filelist_null0=list()
                for(dirpath2,dirnames2,filenames2) in walk(inpath2):
                    filelist_null0.extend(filenames2)
                filelist_null=[fn for fn in filelist_null0 if fn[:12]=='parc_connmat']
                
            elif which_parc==1:#case 2
                inpath2='/imaging/rf02/parcellation/extended_simulations/null_sims40_offset0_bothsnrs/aparc_mod/'
                filelist_null0=list()
                for(dirpath2,dirnames2,filenames2) in walk(inpath2):
                    filelist_null0.extend(filenames2)
                filelist_null=[fn for fn in filelist_null0 if fn[:12]=='parc_connmat']
            elif which_parc==2:#case 2
                inpath2='/imaging/rf02/parcellation/extended_simulations/null_sims40_offset0_bothsnrs/aparc2009/'
                filelist_null0=list()
                for(dirpath2,dirnames2,filenames2) in walk(inpath2):
                    filelist_null0.extend(filenames2)
                filelist_null=[fn for fn in filelist_null0 if fn[:12]=='parc_connmat']
            elif which_parc==3:#case 2
                inpath2='/imaging/rf02/parcellation/extended_simulations/null_sims40_offset0_bothsnrs/aparc2009_mod/'
                filelist_null0=list()
                for(dirpath2,dirnames2,filenames2) in walk(inpath2):
                    filelist_null0.extend(filenames2)
                filelist_null=[fn for fn in filelist_null0 if fn[:12]=='parc_connmat']
            elif which_parc==4:#case 2
                inpath2='/imaging/rf02/parcellation/extended_simulations/null_sims40_offset0_bothsnrs/RG/'
                filelist_null0=list()
                for(dirpath2,dirnames2,filenames2) in walk(inpath2):
                    filelist_null0.extend(filenames2)
                filelist_null=[fn for fn in filelist_null0 if fn[:12]=='parc_connmat']
            
            nullconnmat=np.zeros((nnodes,nnodes,8,nsubs,nsims));
            for fcnt,filename_null in enumerate(filelist_null):#=1:nsims:
                filepath=inpath2+filename_null
                conpre=scipy.io.loadmat(filepath,mat_dtype=True); 
                nullconnmat[:,:,:,:,fcnt]=np.squeeze(conpre['parc'+str(which_parc)+'_connmat'])

            nullconnmat=np.mean(np.abs(nullconnmat),axis=4);
#            base_connmat=repmat(nullconnmat(:,:,[2,1,1,4,3,6,5,5,8,7],:),[1,1,1,1,nsims]);

            basecorr_connmat=np.abs(connmat)#-abs(base_connmat);%ideal 3COH,8imCOH            
            for cc1 in range(10):
                for cc2 in range(nsubs):
                    for cc3 in range(nsims):
                        this_connmat=basecorr_connmat[:,:,cc1,cc2,cc3].copy()+np.triu(1000*np.ones((nnodes,nnodes)))
                        tcm=np.median(this_connmat[this_connmat<1000])
                        this_connmat=this_connmat-np.triu(1000*np.ones((nnodes,nnodes)))
                        this_connmat=this_connmat-tcm
                        basecorr_connmat[:,:,cc1,cc2,cc3]=this_connmat.copy()
            pval=np.zeros((nnodes,nnodes,10,nsims));
            print 'ttest'
            for cc1 in range(10):
                print cc1
                
                for cc2 in range(nsims):
                    this_grandmat=np.squeeze(basecorr_connmat[:,:,cc1,cc2,:])
                    ntests=nnodes*(nnodes-1)/2
                    X=np.zeros((nsubs,ntests))
                    for cc3 in range(nsubs):
                        thismat=this_grandmat[:,:,cc3].copy()+np.triu(1000*np.ones((nnodes,nnodes)))
                        X[cc3,:]=thismat[thismat<1000]
                    tval,thispval,H0=mne.stats.permutation_t_test(X,tail=1)
                    thismat[thismat<1000]=thispval.copy()
                    pval[:,:,cc1,cc2]=thismat.copy()
                   
#                    thisvect=np.squeeze(basecorr_connmat[cc1,cc2,cc3,:,:]);
#                    this_tval=mne.stats.ttest_1samp_no_p(thisvect, sigma=1e-3, method='relative')
#                    this_tval[this_tval<0]=1
#                    for ttcnt,thist in enumerate(this_tval):
#                        pval[cc1,cc2,cc3,ttcnt]=scistats.t.sf(np.abs(thist),nsubs-1) 
            out_name=inpath+'perm_pval_parc'+str(which_parc)+'.mat'
            scipy.io.savemat(out_name,{'pval':pval})
                            
tfin=time.time()-t
print tfin
            
#
#        in_path='/imaging/rf02/parcellation/extended_simulations/sims'+str(n_sims)+'_offset'+str(sim_offset)+'_seeds'+str(nseeds)+'_conns'+str(connperc)+'_bothsnrs_ext'+str(ext)+'/missedcons'
#        with open ('out_path', 'rb') as fp:
#    itemlist = pickle.load(fp)
#        outname=['/aparc/','/aparc_mod/','/aparc2009/','/aparc2009_mod/','/RG/']
#        for thisname in outname:
#            if not os.path.exists(out_path+thisname):
#                os.makedirs(out_path+thisname)
#        
#        this_pathto=out_path+'/activesources_grandmat_sim'+str(perm_inds[0]+sim_offset)+'.py'#'DKA68_orig_modified_coh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot            
#        with open(this_pathto,'wb') as fp:
#            pickle.dump(actsrc_grandmat,fp)
#        this_pathto=out_path+'/activesources_howmanyovs_sim'+str(perm_inds[0]+sim_offset)+'.py'#'DKA68_orig_modified_coh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot            
#        with open(this_pathto,'wb') as fp:
#            pickle.dump(actsrc_ovs,fp)        
#        this_pathto=out_path+'/activesourcesfinal_grandmat_sim'+str(perm_inds[0]+sim_offset)+'.py'#'DKA68_orig_modified_coh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot            
#        with open(this_pathto,'wb') as fp:
#            pickle.dump(actsrcfinal_grandmat,fp)
#        this_pathto=out_path+'/missedcons_grandmat_sim'+str(perm_inds[0]+sim_offset)+'.py'#'DKA68_orig_modified_coh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot            
#        with open(this_pathto,'wb') as fp:
#            pickle.dump(missedcons_grandmat,fp)
#        this_pathto=out_path+'/misscatch_grandmat_sim'+str(perm_inds[0]+sim_offset)+'.py'#'DKA68_orig_modified_coh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot            
#        with open(this_pathto,'wb') as fp:
#            pickle.dump(misscatch_grandmat,fp)
    
#with open ('outfile', 'rb') as fp:
#    itemlist = pickle.load(fp)
# First, we reorder the labels based on their location in the left hemi
#tfin=time.time()-t0
#label_names = [label.name for label in labels]
#
#lh_labels = [name for name in label_names if name.endswith('lh')]
#
## Get the y-location of the label
#label_ypos = list()
#for name in lh_labels:
#    idx = label_names.index(name)
#    ypos = np.mean(labels[idx].pos[:, 1])
#    label_ypos.append(ypos)
#
## Reorder the labels based on their location
#lh_labels = [label for (yp, label) in sorted(zip(label_ypos, lh_labels))]
#
## For the right hemi
#rh_labels = [label[:-2] + 'rh' for label in lh_labels]
#
## Save the plot order and create a circular layout
#node_order = list()
#node_order.extend(lh_labels[::-1])  # reverse the order
#node_order.extend(rh_labels)
#
#node_angles = circular_layout(label_names, node_order, start_pos=90,
#                              group_boundaries=[0, len(label_names) / 2])
#
## Plot the graph using node colors from the FreeSurfer parcellation. We only
## show the 300 strongest connections.
#plt.close("all")
#con_res_mean=np.squeeze(np.mean(con_res_final,axis=2))
#plot_connectivity_circle(np.abs(con_res_mean), range(len(label_names)),n_lines=None,
#                         node_angles=node_angles, node_colors=label_colors,#vmin=0,vmax=1,
#                         title='All-to-All Connectivity '
#                               'Coh')
##con_mod_mean=np.mean(con_mod_final,axis=2)                             
##plot_connectivity_circle(np.abs(con_mod_mean), range(len(label_names)),n_lines=300,
##                         node_angles=node_angles, node_colors=label_colors,#vmin=0,vmax=1e-19,
##                         title='All-to-All Connectivity '
##                               'Coh')
#con_res_mean_ideal=np.squeeze(np.mean(con_res_final_ideal,axis=2))
#plot_connectivity_circle(np.abs(con_res_mean_ideal), label_names, n_lines=None,#vmin=0,#vmax=1,
#                         node_angles=node_angles, node_colors=label_colors,
#                         title='All-to-All Connectivity no leakage Coh')
####
#con_res_mean_realideal=np.squeeze(np.mean(con_res_final_realideal,axis=2))
#plot_connectivity_circle(np.abs(con_res_mean_realideal), label_names, n_lines=None,#vmin=0,#vmax=1,
#                         node_angles=node_angles, node_colors=label_colors,
#                         title='All-to-All Connectivity '
#                               'Ideal Coh')                               
####
#con_res_mean2=np.squeeze(np.mean(con_res_final2,axis=2))
#plot_connectivity_circle(np.abs(con_res_mean2), label_names, n_lines=None,
#                         node_angles=node_angles, node_colors=label_colors,#vmin=0,vmax=1,
#                         title='All-to-All Connectivity '
#                               'ImCoh')
##con2_mod_mean=np.mean(con2_mod_final,axis=2)                             
##plot_connectivity_circle(np.abs(con2_mod_mean), range(len(label_names)),n_lines=300,
##                         node_angles=node_angles, node_colors=label_colors,#vmin=0,vmax=1,
##                         title='All-to-All Connectivity '
##                               'Coh')
#con_res_mean_ideal2=np.squeeze(np.mean(con_res_final_ideal2,axis=2))
#plot_connectivity_circle(np.abs(con_res_mean_ideal2), label_names, n_lines=None,#vmin=0,vmax=1,
#                         node_angles=node_angles, node_colors=label_colors,
#                         title='All-to-All Connectivity '
#                               'no leakage ImCoh') 
#                               
####
#con_res_mean_realideal2=np.squeeze(np.mean(con_res_final_realideal2,axis=2))
#plot_connectivity_circle(np.abs(con_res_mean_realideal2), label_names, n_lines=None,#vmin=0,vmax=1,
#                         node_angles=node_angles, node_colors=label_colors,
#                         title='All-to-All Connectivity '
#                               'Ideal ImCoh')                                
##plt.savefig('circle.png', facecolor='black')
#outname=['aparc','aparc_mod','aparc2009','aparc2009_mod','RG']
##out_path='/imaging/rf02/parcellation/extended_sim/'
##path_to=label_path+outname[which_parc]+'_coh_test.mat'#'DKA68_orig_modified_coh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot
##scipy.io.savemat(path_to,{'con_res_final':con_res_final})
###path_to=label_path+outname[which_parc]+'_coh_mod_test.mat'#'DKA68_orig_modified_coh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot
###scipy.io.savemat(path_to,{'con_mod_final':con_mod_final})
##
##path_to=label_path+outname[which_parc]+'_ideal_coh_test.mat'#'DKA68_ideal_modified_coh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot
##scipy.io.savemat(path_to,{'con_res_final_ideal':con_res_final_ideal})
##
##path_to=label_path+outname[which_parc]+'_imcoh_test.mat'#'DKA68_orig_modified_imcoh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot
##scipy.io.savemat(path_to,{'con_res_final2':con_res_final2})
###path_to=label_path+outname[which_parc]+'_imcoh_mod_test.mat'#'DKA68_orig_modified_imcoh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot
###scipy.io.savemat(path_to,{'con2_mod_final':con2_mod_final})
##
##path_to=label_path+outname[which_parc]+'_ideal_imcoh_test.mat'#'DKA68_ideal_modified_imcoh.mat'#lh.fsaverage_finalpostparc_smooth_split_merge.annot
##scipy.io.savemat(path_to,{'con_res_final_ideal2':con_res_final_ideal2})
#plt.show()
# 
#def mean2(x):
#    y=np.sum(x)/np.size(x)
#    return y
#def corr2(a,b):
#    a=a-mean2(a)
#    b=b-mean2(b)
#    r=(a*b).sum()/np.sqrt((a*a).sum()*(b*b).sum())
#    return r       
#cc2=(con_res_mean_realideal+con_res_mean_realideal.T)/2
##cc2=cc2+np.eye(cc2.shape[0])
#cc1=(con_res_mean+con_res_mean.T)/2 
##cc1=cc1+np.eye(cc1.shape[0])  
#R=corr2(cc1,cc2)