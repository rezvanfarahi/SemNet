# -*- coding: utf-8 -*-
"""
Created on Sat Dec 19 22:44:56 2015

@author: rf02
"""

#    label_names_pre=[label.name for label in labels]
#    labelss_semifinal = [labels_semifinal[jj].morph(subject_from=subject_from, subject_to='fsaverage', subjects_dir=data_path) for jj in range(len(labels_semifinal))]
#    labelss_semifinal_smooth=[labelss_semifinal[jj].smooth(subject='fsaverage', smooth=2, grade=5, subjects_dir=data_path, copy=True) for jj in range(len(labelss_semifinal))]
#    mne.write_labels_to_annot(labelss_semifinal_smooth, subject=subject_from, parc=subject_from+'_parc_morphed_smooth', subjects_dir=data_path,colormap='hsv',hemi='both',overwrite=True)
#    
#    label1v=labels[label_names_pre.index('posterior_inferiorfrontal_caudalmiddlefrontal-lh')].vertices
#    label2v=labels[label_names_pre.index('supramarginal-lh')].vertices
#    vnov=np.in1d(label2v,label1v).sum()
#    print vnov
#    label2ver_new=label2v[vnov]
#     for each individual
#    bb=stc_ctf_mne.in_label(labels)
#    cc=np.abs(stc_ctf_mne.data[:,0])
#    dd=cc.argsort()[::-1]
#    bbv=bb.vertices[0]
#    bbv=np.searchsorted(stc_ctf_mne.vertno[0], bbv)
#    ee=dd[:bb.shape[0]]
#    ff=np.intersect1d(ee,bbv)
#    stc_to.save(data_out)
#     for average
#    stc_to=stc_ctf_mne
#    
#    ###############
#    Matx=stc_to.data
#    stc_all[cnt]=stc_to
#    #Mat_avg=Matx.mean(axis=2)
#    ### combine labels
#    overlap_sub=np.zeros((len(labelss),len(labelss)))
#    for kk, lblk in enumerate(labelss):
#        print kk
#        cc=np.abs(Matx[:,kk])
#        best_vert=sorted(cc,reverse=True)[0]
#        best_verts=np.where(cc>best_vert/2.)[0]
#        
#        for kk2, lblk2 in enumerate(labelss):
#            bb=stc_to.in_label(labelss[kk2])        
#            if kk2%2==0:
#                bbv=bb.vertices[0]
#            else:
#                bbv=bb.vertices[1]+len(stc_to.lh_vertno)
#            
#            overlap_sub[kk,kk2]=len(np.intersect1d(best_verts,bbv))
#    overlap_sub_winner[:,cnt]=np.argmax(overlap_sub, axis=1)
#            
#    print rat[:,cnt].mean()
#from scipy.stats import mode
#vertex_winner=mode(assign_vertex, axis=1)
##label_winner=mode(overlap_sub_winner, axis=1)
#path_c=data_path+'aparc_24sep15_vertex_winner.mat'
##path_c2=data_path+'aparc_14sep15_label_winner.mat'
#
#scipy.io.savemat(path_c,{'vertex_winner':vertex_winner})
##scipy.io.savemat(path_c2,{'label_winner':label_winner})
#
#matx_stc = mne.SourceEstimate(vertex_winner[0], vertices=vertices_to,tmin=1e-3 * 1, tstep=1e-3 * 1, subject='fsaverage')
#path3=data_path+'24sep15_ctf_vertex_winner'
#matx_stc.save(path3)
#unver=unique(assign_vertex[:,0])
#labvernum=np.zeros(unver.shape[0])
#for uc in range(len(unver)):
#    labvernum[uc]=np.where(assign_vertex[:,0]==unver[uc])[0].shape[0]
    
#rat_new=rat.copy()
#rat_blnc[::2,:]=np.maximum(rat[::2,:],rat[1::2,:])
#rat_blnc[1::2,:]=np.maximum(rat[::2,:],rat[1::2,:])
#rat_new=rat_blnc.mean(axis=1)
#labl_new=[labl for ii,labl in enumerate(labelss) if rat_new[ii]>=rat_new.mean()-rat_new.std()/2.]#rat_new.mean()]#-rat_new.std()] 
#annot_fname=data_path+'fsaverage/label/rezvan_parc_cnew_mne_ctf_mn.annot'   
#mne.label.write_labels_to_annot(labl_new, subject='fsaverage', parc='rezvan_parc_cnew_mne_ctf_mn', overwrite=True,
#                          subjects_dir=data_path, hemi='both', verbose=None)
#path_rat=data_path+'icaanalysis_results/mne_ratio_cnew2_mat.mat'
#scipy.io.savemat(path_rat,{'rat':rat})
#Matxl=Matxl.transpose()
#Matxr=Matxr.transpose()
#t_Matxl= mne.stats.ttest_1samp_no_p(Matxl, sigma=0, method='relative')
#t_Matxr= mne.stats.ttest_1samp_no_p(Matxr, sigma=0, method='relative')
#t_Matx=np.zeros((20484,2))
#t_Matx[:,0]=t_Matxl
#t_Matx[:,1]=t_Matxr
#fsave_vertices = [np.arange(10242), np.arange(10242)]
#t_stc = mne.SourceEstimate(t_Matx, vertices=fsave_vertices, tmin=1e-3*2, tstep=1e-3*2, subject='fsaverage', verbose=None)
#out_filet=data_path + 'UVttest_CTF_' + lname 
#t_stc.save(out_filet)
#Ml=Matxl.mean(axis=0)
#Mr=Matxr.mean(axis=0)
#M=np.vstack((Ml,Mr)).transpose()
#t_stc = mne.SourceEstimate(M, vertices=fsave_vertices, tmin=1e-3*2, tstep=1e-3*2, subject='fsaverage', verbose=None)
#out_filet=data_path + 'GrandAverage_CTF_' + lname 
#t_stc.save(out_filet)
#	# stc_psf_meg.save('psf_meg')
#"""
#	from mayavi import mlab
#	fmin = 0.
#	time_label = "EEGMEG %d"
#	fmax = stc_psf_eegmeg.data[:, 0].max()
#	fmid = fmax / 2.
#	brain_eegmeg = stc_psf_eegmeg.plot(surface='inflated', hemi='rh',
#		                           subjects_dir=subjects_dir,
#		                           time_label=time_label, fmin=fmin,
#		                           fmid=fmid, fmax=fmax,
#		                           figure=mlab.figure(size=(500, 500)))
#
#	time_label = "MEG %d"
#	fmax = stc_psf_meg.data[:, 0].max()
#	fmid = fmax / 2.
#	brain_meg = stc_psf_meg.plot(surface='inflated', hemi='rh',
#		                     subjects_dir=subjects_dir,
#		                     time_label=time_label, fmin=fmin,
#		                     fmid=fmid, fmax=fmax,
#		                     figure=mlab.figure(size=(500, 500)))
#
#	# The PSF is centred around the right auditory cortex label,
#	# but clearly extends beyond it.
#	# It also contains "sidelobes" or "ghost sources"
#	# in middle/superior temporal lobe.
#	# For the Aud-RH example, MEG and EEGMEG do not seem to differ a lot,
#	# but the addition of EEG still decreases point-spread to distant areas
#	# (e.g. to ATL and IFG).
#	# The chosen labels are quite far apart from each other, so their PSFs
#	# do not overlap (check in mne_analyze)
#"""
