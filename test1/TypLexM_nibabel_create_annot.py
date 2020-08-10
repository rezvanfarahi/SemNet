# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 23:12:51 2015

@author: rf02
"""
#import sys
#sys.path.append('/imaging/rf02/nibabel')
#
#import nibabel
#import nibabel.freesurfer.io as nifsio
#[aa,bb,cc]=nifsio.read_annot('/imaging/rf02/TypLexMEG/avgsubj/label/lh.aparc.annot')
#atl=nifsio.read_label('/imaging/rf02/TypLexMEG/icaanalysis_results/labels/ATL2-lh.label')
#aa[atl]=2
#cc[2]='ATL'
#bb[2,4]=len(atl)
#nifsio.write_annot('/imaging/rf02/TypLexMEG/avgsubj/labels/lh.mine.annot',aa,bb,cc)

# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 23:12:51 2015

@author: rf02
"""

#mri_annotation2label --subject fsaverage --annotation aparc --hemi lh --outdir /imaging/rf02/TypLexMEG/fsaverage/label
import sys
sys.path.remove('/imaging/local/software/python_packages/nibabel/1.3.0')
sys.path.append('/imaging/rf02')
sys.path.insert(1,'/imaging/local/software/mne_python/latest_v0.9')
import mne
from mne.label import _read_annot as mnera
from mne.label import _write_annot as mnewa

data_path='/imaging/rf02/TypLexMEG'
#import os
#os.chdir('/imaging/rf02')
import nibabel.freesurfer.io as nifsio
[aal,bbl,ccl]=nifsio.read_annot('/imaging/rf02/TypLexMEG/fsaverage/label/lh.aparc.annot')
[aar,bbr,ccr]=nifsio.read_annot('/imaging/rf02/TypLexMEG/fsaverage/label/rh.aparc.annot')
aal,bbl,ccl=mnera('/imaging/rf02/TypLexMEG/fsaverage/label/lh.aparc.annot')
aar,bbr,ccr=mnera('/imaging/rf02/TypLexMEG/fsaverage/label/rh.aparc.annot')



#atl=nifsio.read_label('/imaging/rf02/TypLexMEG/icaanalysis_results/labels/ATL2-lh.label')
#aa[atl]=2
#cc[2]='ATL'
#bb[2,4]=len(atl)

atlr=mne.read_label('/imaging/rf02/TypLexMEG/fsaverage/label/rightATL-rh.label',subject='fsaverage')
atll=mne.read_label('/imaging/rf02/TypLexMEG/fsaverage/label/leftATL-lh.label',subject='fsaverage')
atll_v=atll.vertices#get_vertices_used()
atlr_v=atlr.vertices#get_vertices_used()
aal[atll_v]=2
ccl[2]='ATL'

aar[atlr_v]=2
ccr[2]='ATL'

nifsio.write_annot('/imaging/rf02/TypLexMEG/fsaverage/label/lh.mine.aparc.annot',aal,bbl,ccl)
nifsio.write_annot('/imaging/rf02/TypLexMEG/fsaverage/label/rh.mine.aparc.annot',aar,bbr,ccr)

sfl=mne.read_label('/imaging/rf02/TypLexMEG/fsaverage/label/lh.superiorfrontal.label',subject='fsaverage')
sfl_parts=sfl.split(parts=('front','back'), subject='fsaverage', subjects_dir=data_path)
mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/lh.back_superiorfrontal.label', sfl_parts[0])
mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/lh.front_superiorfrontal.label', sfl_parts[1])

sfr=mne.read_label('/imaging/rf02/TypLexMEG/fsaverage/label/rh.superiorfrontal.label',subject='fsaverage')
sfr_parts=sfr.split(parts=('front','back'), subject='fsaverage', subjects_dir=data_path)
mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/rh.back_superiorfrontal.label', sfr_parts[0])
mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/rh.front_superiorfrontal.label', sfr_parts[1])

sfr_parts=sfr.split(parts=('posterior','middle', 'anterior'), subject='fsaverage', subjects_dir=data_path)

mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/rh.posterior_superiorfrontal.label', sfr_parts[0])
mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/rh.middle_superiorfrontal.label', sfr_parts[1])
mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/rh.anterior_superiorfrontal.label', sfr_parts[2])

sfr=mne.read_label('/imaging/rf02/TypLexMEG/fsaverage/label/rh.precentral.label',subject='fsaverage')
sfr_parts=sfr.split(parts=('superior','middle', 'inferior'), subject='fsaverage', subjects_dir=data_path)
mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/rh.superior_precentral.label', sfr_parts[0])
mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/rh.middle_precentral.label', sfr_parts[1])
mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/rh.inferior_precentral.label', sfr_parts[2])

sfl=mne.read_label('/imaging/rf02/TypLexMEG/fsaverage/label/lh.precentral.label',subject='fsaverage')
sfl_parts=sfl.split(parts=('superior','middle', 'inferior'), subject='fsaverage', subjects_dir=data_path)
mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/lh.superior_precentral.label', sfl_parts[0])
mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/lh.middle_precentral.label', sfl_parts[1])
mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/lh.inferior_precentral.label', sfl_parts[2])

sfr=mne.read_label('/imaging/rf02/TypLexMEG/fsaverage/label/rh.postcentral.label',subject='fsaverage')
sfr_parts=sfr.split(parts=('superior','middle', 'inferior'), subject='fsaverage', subjects_dir=data_path)
mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/rh.superior_postcentral.label', sfr_parts[0])
mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/rh.middle_postcentral.label', sfr_parts[1])
mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/rh.inferior_postcentral.label', sfr_parts[2])


sfl=mne.read_label('/imaging/rf02/TypLexMEG/fsaverage/label/lh.postcentral.label',subject='fsaverage')
sfl_parts=sfl.split(parts=('superior','middle', 'inferior'), subject='fsaverage', subjects_dir=data_path)
mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/lh.superior_postcentral.label', sfl_parts[0])
mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/lh.middle_postcentral.label', sfl_parts[1])
mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/lh.inferior_postcentral.label', sfl_parts[2])


sfr=mne.read_label('/imaging/rf02/TypLexMEG/fsaverage/label/rh.fusiform.label',subject='fsaverage')
sfr_parts=sfr.split(parts=('anterior','posterior'), subject='fsaverage', subjects_dir=data_path)
mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/rh.posterior_fusiform.label', sfr_parts[0])
mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/rh.anterior_fusiform.label', sfr_parts[1])


sfl=mne.read_label('/imaging/rf02/TypLexMEG/fsaverage/label/lh.fusiform.label',subject='fsaverage')
sfl_parts=sfl.split(parts=('anterior','posterior'), subject='fsaverage', subjects_dir=data_path)
mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/lh.posterior_fusiform.label', sfl_parts[0])
mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/lh.anterior_fusiform.label', sfl_parts[1])

### NOTE!! I did this with my modified label split version to take the second eigen vector sys.path.insert(1,'/imaging/rf02/TypLexMEG/mne_python0.9')
sfl=mne.read_label('/imaging/rf02/TypLexMEG/fsaverage/label/lh.bankssts.label',subject='fsaverage')
sfl_parts=sfl.split(parts=('superior', 'inferior'), subject='fsaverage', subjects_dir=data_path)
mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/lh.superior_bankssts.label', sfl_parts[1])
mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/lh.inferior_bankssts.label', sfl_parts[0])

sfr=mne.read_label('/imaging/rf02/TypLexMEG/fsaverage/label/rh.bankssts.label',subject='fsaverage')
sfr_parts=sfr.split(parts=('superior', 'inferior'), subject='fsaverage', subjects_dir=data_path)
mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/rh.superior_bankssts.label', sfr_parts[1])
mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/rh.inferior_bankssts.label', sfr_parts[0])

labelname_lsb='/imaging/rf02/TypLexMEG/fsaverage/label/lh.superior_bankssts.label'    
labelname_lib='/imaging/rf02/TypLexMEG/fsaverage/label/lh.inferior_bankssts.label'
labelname_lst='/imaging/rf02/TypLexMEG/fsaverage/label/lh.superiortemporal.label'    
labelname_lmt='/imaging/rf02/TypLexMEG/fsaverage/label/lh.middletemporal.label' 
labelname_ltt='/imaging/rf02/TypLexMEG/fsaverage/label/lh.transversetemporal.label'       

label_lsb=mne.read_label(labelname_lsb,subject='fsaverage')
label_lib=mne.read_label(labelname_lib,subject='fsaverage')
label_lst=mne.read_label(labelname_lst,subject='fsaverage')
label_lmt=mne.read_label(labelname_lmt,subject='fsaverage')
label_ltt=mne.read_label(labelname_ltt,subject='fsaverage')


label_lmidtemp=label_lmt+label_lib
mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/lh.middle_temporal.label', label_lmidtemp)
label_lsuptemp=label_lst+label_lsb+label_ltt
mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/lh.superior_temporal.label', label_lsuptemp)


labelname_rsb='/imaging/rf02/TypLexMEG/fsaverage/label/rh.superior_bankssts.label'    
labelname_rib='/imaging/rf02/TypLexMEG/fsaverage/label/rh.inferior_bankssts.label'
labelname_rst='/imaging/rf02/TypLexMEG/fsaverage/label/rh.superiortemporal.label'    
labelname_rmt='/imaging/rf02/TypLexMEG/fsaverage/label/rh.middletemporal.label' 
labelname_rtt='/imaging/rf02/TypLexMEG/fsaverage/label/rh.transversetemporal.label'       

label_rsb=mne.read_label(labelname_rsb,subject='fsaverage')
label_rib=mne.read_label(labelname_rib,subject='fsaverage')
label_rst=mne.read_label(labelname_rst,subject='fsaverage')
label_rmt=mne.read_label(labelname_rmt,subject='fsaverage')
label_rtt=mne.read_label(labelname_rtt,subject='fsaverage')


label_rmidtemp=label_rmt+label_rib
mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/rh.middle_temporal.label', label_rmidtemp)
label_rsuptemp=label_rst+label_rsb+label_rtt
mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/rh.superior_temporal.label', label_rsuptemp)

####ATL
##left
atl=mne.read_label('/imaging/rf02/TypLexMEG/fsaverage/label/lh.parsopercularis.label',subject='fsaverage')
tpl=mne.read_label('/imaging/rf02/TypLexMEG/fsaverage/label/lh.parstriangularis.label',subject='fsaverage')
posterior_inferiorfrontal=atl+tpl
mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/lh.posterior_inferiorfrontal.label', posterior_inferiorfrontal)

##right
atr=mne.read_label('/imaging/rf02/TypLexMEG/fsaverage/label/rh.parsopercularis.label',subject='fsaverage')
tpr=mne.read_label('/imaging/rf02/TypLexMEG/fsaverage/label/rh.parstriangularis.label',subject='fsaverage')
posterior_inferiorfrontal=atr+tpr
mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/rh.posterior_inferiorfrontal.label', posterior_inferiorfrontal)


##left
atl=mne.read_label('/imaging/rf02/TypLexMEG/fsaverage/label/lh.lateralorbitofrontal.label',subject='fsaverage')
tpl=mne.read_label('/imaging/rf02/TypLexMEG/fsaverage/label/lh.parsorbitalis.label',subject='fsaverage')
anterior_inferiorfrontal=atl+tpl
mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/lh.anterior_inferiorfrontal.label', anterior_inferiorfrontal)

##right
atr=mne.read_label('/imaging/rf02/TypLexMEG/fsaverage/label/rh.lateralorbitofrontal.label',subject='fsaverage')
tpr=mne.read_label('/imaging/rf02/TypLexMEG/fsaverage/label/rh.parsorbitalis.label',subject='fsaverage')
anterior_inferiorfrontal=atr+tpr
mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/rh.anterior_inferiorfrontal.label', anterior_inferiorfrontal)


##left

#satl=mne.read_label('/imaging/rf02/TypLexMEG/fsaverage/label/supATL-lh.label',subject='fsaverage')
atl=mne.read_label('/imaging/rf02/TypLexMEG/fsaverage/label/lh.ATL.label',subject='fsaverage')
tpl=mne.read_label('/imaging/rf02/TypLexMEG/fsaverage/label/lh.temporalpole.label',subject='fsaverage')
ATLnew=atl+tpl
mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/lh.ATLnew.label', ATLnew)

##right
#satr=mne.read_label('/imaging/rf02/TypLexMEG/fsaverage/label/supATL-rh.label',subject='fsaverage')
atr=mne.read_label('/imaging/rf02/TypLexMEG/fsaverage/label/rh.ATL.label',subject='fsaverage')
tpr=mne.read_label('/imaging/rf02/TypLexMEG/fsaverage/label/rh.temporalpole.label',subject='fsaverage')
ATLnew=atr+tpr
mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/rh.ATLnew.label', ATLnew)



sfr=mne.read_label('/imaging/rf02/TypLexMEG/fsaverage/label/rh.precentral.label',subject='fsaverage')
sfr_parts=sfr.split(parts=('superior2', 'inferior2'), subject='fsaverage', subjects_dir=data_path)
mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/rh.superior2_precentral.label', sfr_parts[0])
mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/rh.inferior2_precentral.label', sfr_parts[1])

sfl=mne.read_label('/imaging/rf02/TypLexMEG/fsaverage/label/lh.precentral.label',subject='fsaverage')
sfl_parts=sfl.split(parts=('superior2', 'inferior2'), subject='fsaverage', subjects_dir=data_path)
mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/lh.superior2_precentral.label', sfl_parts[0])
mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/lh.inferior2_precentral.label', sfl_parts[1])

sfr=mne.read_label('/imaging/rf02/TypLexMEG/fsaverage/label/rh.postcentral.label',subject='fsaverage')
sfr_parts=sfr.split(parts=('superior2', 'inferior2'), subject='fsaverage', subjects_dir=data_path)
mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/rh.superior2_postcentral.label', sfr_parts[0])
mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/rh.inferior2_postcentral.label', sfr_parts[1])


sfl=mne.read_label('/imaging/rf02/TypLexMEG/fsaverage/label/lh.postcentral.label',subject='fsaverage')
sfl_parts=sfl.split(parts=('superior2', 'inferior2'), subject='fsaverage', subjects_dir=data_path)
mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/lh.superior2_postcentral.label', sfl_parts[0])
mne.write_label('/imaging/rf02/TypLexMEG/fsaverage/label/lh.inferior2_postcentral.label', sfl_parts[1])
