# -*- coding: utf-8 -*-
"""
Created on Mon May  4 23:26:35 2015

@author: rf02
"""
import sys
sys.path.append('/imaging/local/software/mne_python/latest')
import mne
import numpy as np
from mne.label import _read_annot as mnera
from mne.label import _write_annot as mnewa
from mne.label import write_labels_to_annot 
from mne.label import _n_colors as nc
data_path='/imaging/rf02/TypLexMEG'

#labellist=['middle_superiorfrontal','inferior_precentral', 'middle_temporal', 'anterior_inferiorfrontal','anterior_fusiform','superior_precentral', 'inferior_postcentral',  'superior_postcentral', 'anterior_superiorfrontal', 'middle_precentral', 'posterior_fusiform', 'superior_temporal', 'posterior_inferiorfrontal','posterior_superiorfrontal', 'middle_postcentral','ATLnew']
labellist=['middle_superiorfrontal','inferior2_precentral', 'middle_temporal', 'anterior_fusiform','superior2_precentral', 'inferior2_postcentral',  'superior2_postcentral', 'anterior_superiorfrontal',  'posterior_fusiform', 'superior_temporal', 'posterior_inferiorfrontal','posterior_superiorfrontal', 'ATL_latmed', 'AG_new', 'vWFA2']


#annotl,ctabl,namesl=mnera('/imaging/rf02/TypLexMEG/SS1_Meg10_0378/label/lh.SS1_Meg10_0378_parc.annot')#('/imaging/rf02/TypLexMEG/fsaverage/label/lh.aparc.annot')
annotl,ctabl,namesl=mnera('/imaging/rf02/TypLexMEG/fsaverage/label/lh.aparc.annot')
annotr,ctabr,namesr=mnera('/imaging/rf02/TypLexMEG/fsaverage/label/rh.aparc.annot')

CTAB_L=np.zeros((ctabl.shape[0]+len(labellist),ctabl.shape[1]))  
CTAB_R=np.zeros((ctabr.shape[0]+len(labellist),ctabr.shape[1])) 

colors=nc(ctabl.shape[0]+len(labellist),bytes_=True, cmap='gist_ncar')
annot_id_coding = np.array((1, 2 ** 8, 2 ** 16))
annot_id=np.sum(colors[:,:3] * annot_id_coding, axis=1)
CTAB_L[:,:4]=colors; CTAB_L[:,4]=annot_id
CTAB_R[:,:4]=colors; CTAB_R[:,4]=annot_id
CTAB_L[:ctabl.shape[0],:]=ctabl
CTAB_R[:ctabr.shape[0],:]=ctabr

NAMES_L=namesl+labellist
NAMES_R=namesr+labellist

for ii, label in enumerate(labellist):
    print label
    labelname_l='/imaging/rf02/TypLexMEG/fsaverage/label/lh.' + label + '.label'    
    labelname_r='/imaging/rf02/TypLexMEG/fsaverage/label/rh.' + label + '.label'    

    label_l=mne.read_label(labelname_l,subject='fsaverage')
    label_r=mne.read_label(labelname_r,subject='fsaverage')
    annotl[label_l.vertices]=CTAB_L[-len(labellist)+ii,4]
    annotr[label_r.vertices]=CTAB_R[-len(labellist)+ii,4]
    
['bankssts', 'corpuscallosum', 'fusiform', 'lateralorbitofrontal',  'parsopercularis', 'parsorbitalis', 'parstriangularis',  'postcentral',  'precentral',  'superiorfrontal', 'superiortemporal', 'supramarginal', 'frontalpole', 'temporalpole', 'transversetemporal', 'insula', 'anterior_inferiorfrontal', 'posterior_inferiorfrontal', 'inferior_precentral', 'middle_precentral', 'superior_precentral', 'inferior_postcentral', 'middle_postcentral', 'superior_postcentral', 'anterior_superiorfrontal', 'middle_superiorfrontal', 'posterior_superiorfrontal', 'anterior_fusiform', 'posterior_fusiform', 'superior_temporal', 'middle_temporal', 'ATLsmall']

#write_labels_to_annot(labels_l, subject=None, parc=None, overwrite=False, subjects_dir=None, annot_fname='lh.rez_aparc', colormap='hsv', hemi='lh', verbose=None)
fname='/imaging/rf02/TypLexMEG/fsaverage/label/lh.rez_aparc_25oct16.annot'
mnewa(fname, annotl, CTAB_L, NAMES_L)

fname='/imaging/rf02/TypLexMEG/fsaverage/label/rh.rez_aparc_25oct16.annot'
mnewa(fname, annotr, CTAB_R, NAMES_R)
