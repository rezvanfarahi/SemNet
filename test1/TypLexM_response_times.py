# -*- coding: utf-8 -*-
"""
Created on Fri May 15 11:58:50 2015

@author: rf02
"""
1
1066.65
	

939.8452
	

1003.248

2
	

850.2529
	

705.6374
	

777.9451

3
	

940.5432
	

738.8222
	

839.6827

4
	

826.7209
	

737.8901
	

782.3055

5
	

879.6824
	

804.4091
	

842.0457
			
			

8!!!6
	

873.2976
	

842.1125
	

857.7051

9
	

804.3977
	

681.1739
	

742.7858

10
	

972.4176
	

909.5604
	

940.989
			

12
	

617.8193
	

554.8256
	

586.3224

13
	

921.4944
	

808.5114
	

865.0029

14
	

673.7011
	

596.3556
	

635.0284

15
	

934.9011
	

787.0111
	

860.9561

16
	

933.8718
	

749.7253
	

841.7985
			

18
	

885.4023
	

872.5595
	

878.9809

19
	

1018.488
	

794.8222
	

906.6552

20
	

752.9535
	

742.5455
	

747.7495

21
	

1003.296
	

971.6145
	

987.4554
rt_c=np.array([939.8452,705.6374, 738.8222, 737.8901, 804.4091,842.1125, 681.1739, 909.5604,554.8256,808.5114,596.3556, 787.0111, 749.7253, 872.5595, 794.8222, 742.5455,971.6145 ])
data_out=data_path+'RT_concrete'
scipy.io.savemat(data_out,{'rt_c':rt_c})
rt_a=np.array([1066.65, 850.2529, 940.5432, 826.7209, 879.6824, 873.2976, 804.3977, 972.4176, 617.8193, 921.4944, 673.7011,934.9011,933.8718, 885.4023, 1018.488,752.9535, 1003.296])
data_out=data_path+'RT_Abstract'
scipy.io.savemat(data_out,{'rt_a':rt_a})
#RT=np.hstack([rt_c,rt_a])#; RT=1/RT
RT=np.subtract(rt_c,rt_a); RT=RT/np.max(np.abs(RT))
condition=np.linspace(0,1,RT.shape[0])
intercept=np.ones((RT.shape[0],),dtype=np.float)
names=['RT','condition']#'intercept',, 'interaction'
#condition=intercept.copy()
#condition[18:]=-1
interaction=RT*condition
design_matrix=np.column_stack([RT,condition])#intercept,,interaction
lm=mne.stats.linear_regression(stcs,design_matrix, names)
lmc=lm['condition']
lmcp=lmc.p_val
lmr=lm['RT']
lmi=lm['intercept']

labelss = mne.read_labels_from_annot(subject='fsaverage', parc='rezvan_annot', subjects_dir=data_path)
labelss.remove(labelss[-2])
labelss.remove(labelss[-1])
bb=lmcp.in_label(np.sum(labelss))
fsave_vertices = [np.arange(10242), np.arange(10242)]
nnl=np.in1d(fsave_vertices[0],bb.lh_vertno)
nnr=np.in1d(fsave_vertices[1],bb.rh_vertno)
include_ver=np.hstack((fsave_vertices[0][nnl], fsave_vertices[0][nnr]+10242))
lmcp=lmr.p_val
#MTs=MT.reshape(MT.shape[1]*MT.shape[0],)
pval=lmcp.data#[include_ver,:]
r,p=mne.stats.fdr_correction(pval)
p.min()
pcorr=1-p
tmin1=0
tstep1=50
vertices_to = [np.arange(10242), np.arange(10242)]
pmat=np.zeros(lmcp.data.shape)
pmat[include_ver,:]=pcorr
matx_stc = mne.SourceEstimate(pmat, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')

#datag=np.log(stc_grand.data)
#datag_stc = mne.SourceEstimate(datag, vertices=stc_grand.vertno, tmin=stc_grand.tmin, tstep=stc_grand.tstep, subject='fsaverageect')
matx_stc=lmc.t_val
out_file=out_path + 'Cregression_firstMorphed_SemLoc_icaclean_Subtract_Source_Evoked_m500_700'
matx_stc.save(out_file)

ppp=((np.arange(pval.shape[0])+1)*0.05)/(pval.shape[0])
pval.sort()
fd=np.subtract(pval,ppp)
fdf=np.where(fd<=0)
fdf
pval[fdf[0][0]]