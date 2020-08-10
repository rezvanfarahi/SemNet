# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 23:54:29 2015

@author: rf02
"""
import mne
import numpy as np
import scipy
labelss = mne.read_labels_from_annot(subject='fsaverage', parc='rezvan_annot', subjects_dir=data_path)
labelss.remove(labelss[-2])
labelss.remove(labelss[-1])
bb=lmcp.in_label(np.sum(labelss))
fsave_vertices = [np.arange(10242), np.arange(10242)]
nnl=np.in1d(fsave_vertices[0],bb.lh_vertno)
nnr=np.in1d(fsave_vertices[1],bb.rh_vertno)
include_ver=np.hstack((fsave_vertices[0][nnl], fsave_vertices[0][nnr]+10242))

#MTs=MT.reshape(MT.shape[1]*MT.shape[0],)
MTs=MT[include_ver,9]
pval=scipy.stats.t.sf(np.abs(np.unique(MTs[:,0])),16)*2
ppp=((np.arange(pval.shape[0])+1)*0.05)/(pval.shape[0])
pval.sort()
fd=np.subtract(pval,ppp)
fdf=np.where(fd<=0)
fdf
pval[fdf[0][0]]
