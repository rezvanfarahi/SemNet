# -*- coding: utf-8 -*-
"""
Created on Mon May  4 23:26:35 2015

@author: rf02
"""
import sys
sys.path.insert(1,'/imaging/local/software/mne_python/latest')
import mne
import numpy as np
from mne.label import _read_annot as mnera
from mne.label import _write_annot as mnewa
from mne.label import write_labels_to_annot 
from mne.label import _n_colors as nc
data_path='/imaging/rf02/TypLexMEG'

labellist=['anterior_inferiorfrontal','inferior_precentral', 'posterior_inferiorfrontal','middle_precentral', 'inferior_postcentral', 'middle_postcentral', 'anterior_superiorfrontal','superior_postcentral', 'middle_superiorfrontal','superior_precentral',  'posterior_superiorfrontal', 'anterior_fusiform',  'superior_temporal', 'middle_temporal','posterior_fusiform','ATLnew']

annotl,ctabl,namesl=mnera('/imaging/rf02/TypLexMEG/fsaverage/label/lh.aparc.annot')
annotr,ctabr,namesr=mnera('/imaging/rf02/TypLexMEG/fsaverage/label/rh.aparc.annot')


#CTAB_L[:ctabl.shape[0],:]=ctabl
#CTAB_R[:ctabr.shape[0],:]=ctabr


labelss = mne.read_labels_from_annot(subject='fsaverage', parc='aparc', subjects_dir=data_path)
lnew=2*range(len(labellist))
labelss=labelss[:-1]
namesl=[str(label.name[:-3]) for label in labelss[::2]]
namesr=[str(label.name[:-3]) for label in labelss[1::2]]

NAMES_L=namesl+labellist+['unknown']
NAMES_R=namesr+labellist+['unknown']

CTAB_L=np.concatenate((ctabl,ctabl[:15,:]),axis=0) 
CTAB_R=np.concatenate((ctabr,ctabr[:15,:]),axis=0)

colors=nc(len(NAMES_L),bytes_=True, cmap='Accent')
annot_id_coding = np.array((1, 2 ** 8, 2 ** 16))
annot_id=np.sum(colors[:,:3] * annot_id_coding, axis=1)
CTAB_L[:,:4]=colors; CTAB_L[:,4]=annot_id
CTAB_R[:,:4]=colors; CTAB_R[:,4]=annot_id

for ii, label in enumerate(labellist):
    print label
    labelname_l='/imaging/rf02/TypLexMEG/fsaverage/label/lh.' + label + '.label'    
    labelname_r='/imaging/rf02/TypLexMEG/fsaverage/label/rh.' + label + '.label'    

    label_l=mne.read_label(labelname_l,subject='fsaverage')
    label_r=mne.read_label(labelname_r,subject='fsaverage')
    lnew[2*ii]=label_l
    lnew[2*ii+1]=label_r
labels_all=labelss+lnew
annot_l=annotl-annotl
annot_r=annotr-annotr
for ii in range(len(labels_all)/2):
    annot_l[labels_all[ii*2].vertices]=CTAB_L[ii,4]
    annot_r[labels_all[ii*2+1].vertices]=CTAB_R[ii,4]
annot_l[np.where(annot_l==0)]=CTAB_L[-1,4]
annot_l[np.where(annot_r==0)]=CTAB_R[-1,4]

['bankssts', 'corpuscallosum', 'fusiform', 'lateralorbitofrontal',  'parsopercularis', 'parsorbitalis', 'parstriangularis',  'postcentral',  'precentral',  'superiorfrontal', 'superiortemporal', 'supramarginal', 'frontalpole', 'temporalpole', 'transversetemporal', 'insula', 'anterior_inferiorfrontal', 'posterior_inferiorfrontal', 'inferior_precentral', 'middle_precentral', 'superior_precentral', 'inferior_postcentral', 'middle_postcentral', 'superior_postcentral', 'anterior_superiorfrontal', 'middle_superiorfrontal', 'posterior_superiorfrontal', 'anterior_fusiform', 'posterior_fusiform', 'superior_temporal', 'middle_temporal', 'ATLnew']

#write_labels_to_annot(labels_l, subject=None, parc=None, overwrite=False, subjects_dir=None, annot_fname='lh.rez_aparc', colormap='hsv', hemi='lh', verbose=None)
fname='/imaging/rf02/TypLexMEG/fsaverage/label/lh.rez_aparc_cnew.annot'
mnewa(fname, annot_l, CTAB_L, NAMES_L)

fname='/imaging/rf02/TypLexMEG/fsaverage/label/rh.rez_aparc_cnew.annot'
mnewa(fname, annot_r, CTAB_R, NAMES_R)
