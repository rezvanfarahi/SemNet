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
#sys.path.append('/imaging/local/software/mne_python/latest_v0.9')
sys.path.insert(1,'/imaging/rf02/TypLexMEG/mne_python0.9')
#sys.path.insert(1,'/imaging/local/software/mne_python/latest_v0.9full')

sys.path.insert(1,'/imaging/rf02/scikit-learn-0.16.1')
sys.path.insert(1,'/home/rf02/.local/lib/python2.7/site-packages')


# for qsub
# (add ! here if needed) 
sys.path.append('/imaging/local/software/EPD/latest/x86_64/bin/python')
#sys.path.append('/imaging/rf02/TypLexMEG/mne_python_v8/')
#sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
###
import sklearn
from scipy.stats.mstats import zscore
import os
import numpy as np
import mne
from mne import io
from mne.io import Raw
import operator
from mne.minimum_norm import (read_inverse_operator, make_inverse_operator, write_inverse_operator, apply_inverse, apply_inverse_epochs,psf_ctf_new, point_spread_function,psf_ctf_18Dec15)
from mne.minimum_norm.inverse import (prepare_inverse_operator)
import scipy.io
from mne import find_events
from mne.connectivity import seed_target_indices, spectral_connectivity
from mne import read_labels_from_annot
from mne.label import _n_colors as nc
from mne.label import _read_annot as mnera
from mne.label import _write_annot as mnewa
import copy
import pickle
data_path = '/imaging/rf02/TypLexMEG/' # root directory for your MEG data
os.chdir(data_path)
event_path = '/group/erp/data/olaf.hauk/MEG/TypLexM/data/MF22_new/'
#label_path = '/group/erp/data/olaf.hauk/MEG/TypLexM/labels/createdlabels_SL/'
inv_path = '/imaging/rf02/TypLexMEG/'
subjects_dir = data_path 
print sys.argv
subject_inds = []
for ss in sys.argv[1:]:
   subject_inds.append( int( ss ) )

subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]
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
'fsaverage'
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
'SS21_Meg11_0147',
'fsaverage'
]

ll = []
for ss in subject_inds:
 ll.append(list_all[ss])


print "ll:"
print ll
ii=-1

labelss = mne.read_labels_from_annot(subject='fsaverage', parc='aparc.a2009s', subjects_dir=subjects_dir)#parc='rez_aparc_24sep15'

#labelss.remove(labelss[-3])
labelss.remove(labelss[-1])
labelss.remove(labelss[-1])
print len(labelss)
llk=labelss[0]+labelss[-2]
lrk=labelss[1]+labelss[-1]
llk.name=labelss[0].name
lrk.name=labelss[1].name
#labelss.remove(labelss[-1])
#labelss.remove(labelss[-1])
Matxl=np.zeros((20484,17))
Matxr=np.zeros((20484,17))
Matx=np.zeros((20484,len(labelss)+1,17))
rat=np.zeros((len(labelss),len(labelss),17))
rat_blnc=np.zeros((len(labelss),17))
stc_all=range(17)
overlap_sub_winner=np.zeros((len(labelss),17))
assign_vertex=np.zeros((20484,17))

label_colors = [label.color for label in labelss]
# read label(s)
for lbl in labelss:
    lbl.values.fill(1.0)
for cnt, meg in enumerate(ll):
    ii=ii+1
    print cnt
    subjects_dir = data_path 
    
    vertices_to = [np.arange(10242), np.arange(10242)]
    subject_no=subject_inds[cnt]

    # read inverse operators
    if cnt==17:
        labels = mne.read_labels_from_annot(subject=subjects[subject_no], parc='aparc.a2009s', subjects_dir=subjects_dir)
        labels.remove(labels[-1])
        labels.remove(labels[-1])#only 2009
    else:
        labels = mne.read_labels_from_annot(subject=subjects[subject_no], parc='morphed_aparc_2009_21Dec15', subjects_dir=subjects_dir)#parc='rez_aparc_24sep15'
        fnamem=data_path + subjects[subject_no]+'/label/morphed_aparc_2009_09Jul16'
        with open(fnamem, 'rb') as f:
            labels=pickle.load(f)
    print labels[-1].name
    labelname=[labels[ii1].name for ii1 in range(len(labels))]
    if 'frontalpole-lh' in labelname:
        labels.remove(labels[labelname.index('frontalpole-lh')])
        labels.remove(labels[labelname.index('frontalpole-lh')])#don't panic! it's actually right hemishphere
    #labels = [labelss[jj].morph(subject_from='fsaverage', smooth=5, subject_to=subjects[subject_no], subjects_dir=data_path) for jj in range(len(labelss))]
    labels_forcolor = copy.deepcopy(labels)

    snr = 3.0#ediv.mean()
    lambda2 = 1.0 / snr ** 2
    method = 'MNE'  # can be 'MNE' or 'sLORETA'
    mode = 'svd'
    n_svd_comp = 1
    lblname_cnt=0
    
    
    ## make and write divided annot
    path_c=data_path+'explained90all_2009.mat'
    ltc=scipy.io.loadmat(path_c,mat_dtype=True); explained90=np.squeeze(ltc['EXP90'])
    if 'frontalpole-lh' in labelname:
        explained90=np.delete(explained90,[labelname.index('frontalpole-lh'),labelname.index('frontalpole-lh')+1],0)
#    explained90=np.ones(explained90.shape)
    llnew=[labels[ii2] for ii2 in range(len(labels)) if explained90[ii2]==1 or len(labels[ii2].vertices)<10]#for aparc
#    llnew=[labels[ii2] for ii2 in range(len(labels)) if explained90[ii2]==1] # for aparc2009
    lln1=[labels[ii3].split(explained90[ii3], subjects_dir=data_path)[:] for ii3 in range(len(labels)) if explained90[ii3]>1 and len(labels[ii3].vertices)>=10]
#    lln1=[labels[ii3].split(explained90[ii3], subjects_dir=data_path)[:] for ii3 in range(len(labels)) if explained90[ii3]>1]
    llnn=list()
    for ii in range(len(lln1)):
        llnn+=lln1[ii]
    
    labels_backup=copy.deepcopy(labels)
    labels=llnn+llnew
    label_names=[labels[lcn3].name for lcn3 in range(len(labels))]
    ordered_indices=np.zeros((len(label_names),1))
    lcntl=-1
    for lcnt3 in range(len(label_names)):
        if label_names[lcnt3].endswith('lh'):
            lcntl=lcntl+1
            thislabelnameR=label_names.index(label_names[lcnt3][:-3]+'-rh')
            this_indices=np.array([lcnt3,thislabelnameR])
            ordered_indices[2*lcntl]=lcnt3; ordered_indices[2*lcntl+1]=thislabelnameR
    labels=[labels[jjk] for jjk in np.squeeze(ordered_indices.astype(int))]
    subject_from = subjects[subject_no]
    subject_to = 'fsaverage'
    colors=nc(len(labels),bytes_=False, cmap='gist_ncar')#Note! left right different colors, not important here.
    colors=colors[np.random.permutation(int(len(labels)))]
    for lcnt1 in range(len(labels)):
        labels[lcnt1].color=colors[lcnt1,:]
    fnamemf=data_path + subject_from+'/label/'+subject_from+'_morphed_split_aparc_2009'
    with open(fnamemf,'wb') as f:
        pickle.dump(labels,f)
    if cnt==17:
        mne.write_labels_to_annot(labels, subject=subject_from, parc=subject_from+'_morphed_split_aparc_2009', subjects_dir=data_path,colormap='gist_ncar',hemi='both',overwrite=True)
    