# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 20:05:51 2016

@author: rf02
"""
def remove_fake_R(labels, Rmat, con_mat, psd_mat ):

    #for cnt1 in range(len(labels)): #seed label
    #    for cnt2 in range(len(labels)): #all labels
    #        this_patch=label_patches[cnt2]
    #        for cnt3 in this_patch:
    #            if cnt1!=cnt3:
    #                sub_patch=np.asarray([thisp for thisp in this_patch if thisp!=cnt3])
    #                coef_mat=Rmat[sub_patch,:][:,sub_patch]
    #                Gzz=Cross_spec[cnt1]*np.ones((len(sub_patch),))
    #                scale_coeff0=Gzz*Cross_spec[sub_patch]
    #                scale_coeff=np.divide(coef_mat,scale_coeff0[:,np.newaxis])
    #                con_vect=con_mat[sub_patch,cnt1]
    #                est_corrs=np.linalg.solve(scale_coeff,con_vect)
    #                est_corr=Rmat[cnt3,:]*est_corrs
    #                act_corr=con_mat[cnt1,cnt3]
    #                if (np.abs(est_corr-act_corr)/np.abs(act_corr))<0.1:
    #                    con_mat_modified[cnt1,cnt3]=0
    #                    con_mat_modified[cnt3,cnt1]=0
    # 
#    Rmat=Rmat.transpose() 
    #Rmat: rows leakage to
    for cnt1 in range(len(labels)): #seed label
        for cnt2 in range(len(labels)): #destination labels
            if cnt1<cnt2:
                con_mat[cnt1,cnt2,:]=con_mat[cnt2,cnt1,:]
    
    sRmat=np.divide(Rmat,np.max(Rmat,axis=0))  
    label_patches=np.where(sRmat>0.25)   
             
    for cnt1 in range(len(labels)): #seed label
        for cnt2 in range(len(labels)): #destination labels
            if cnt1!=cnt2:
                this_patch=label_patches[1][label_patches[0]==cnt2]
                thislabw=np.where(this_patch==cnt2)[0]
                coef_mat=Rmat[this_patch,:][:,this_patch].transpose()
                est_corrs=np.zeros(psd_mat.shape[1],1)
                for fcnt in range(psd_mat.shape[1]):
                    Gzz=psd_mat[cnt1,fcnt]*np.ones((len(this_patch),))
                    scale_coeff0=np.sqrt(Gzz*psd_mat[this_patch,fcnt])
                    scale_coeff=np.divide(coef_mat,scale_coeff0[:,np.newaxis])
                    con_vect=con_mat[this_patch,cnt1]
                    est_corrs=np.linalg.solve(scale_coeff,con_vect)
                    est_corr[fcnt]=est_corrs[thislabw]
                con_mat_modified[cnt1,cnt2]=np.mean(est_corr)
                                


def remove_fake_R(labels, Rmat, con_mat, psd_mat ):
    for cnt1 in range(len(labels)): #seed label
        for cnt2 in range(len(labels)): #destination labels
            if cnt1<cnt2:
                con_mat[cnt1,cnt2,:]=con_mat[cnt2,cnt1,:]
    con_mat_mean=np.mean(con_mat,axis=2)
    sRmat=np.divide(Rmat,np.max(Rmat,axis=0))  
    label_patches=np.where(sRmat>0.25)   
    con_mat_modified=np.zeros((len(labels),len(labels)))        
    for cnt1 in range(len(labels)): #seed label
        print cnt1
        for cnt2 in range(len(labels)): #destination labels
            if cnt1!=cnt2:
                this_patch=label_patches[0][label_patches[1]==cnt2]
                thislabw=np.where(this_patch==cnt2)[0]
                notthislabw=np.where(this_patch!=cnt2)[0]
                coef_mat_all=Rmat[this_patch,:][:,this_patch].transpose()
                coef_mat=coef_mat_all[notthislabw,:][:,notthislabw]
                this_coef_vect=coef_mat_all[thislabw,notthislabw]
                est_corrs=np.zeros((psd_mat.shape[1],1))
                for fcnt in range(psd_mat.shape[1]):
                    Gzz=psd_mat[cnt1,fcnt]*np.ones((len(this_patch)-1,))
                    scale_coeff0=np.sqrt(Gzz*psd_mat[this_patch[notthislabw],fcnt])
                    psd_patch=psd_mat[this_patch[notthislabw],fcnt]
                    scale_coeff=np.divide(coef_mat,scale_coeff0[:,np.newaxis])
                    con_vect=con_mat[this_patch[notthislabw],cnt1,fcnt]
                    est_Sxys=np.divide(con_vect,np.diag(scale_coeff))#np.linalg.solve(scale_coeff,con_vect)
                    est_Sxy=np.max(est_Sxys*this_coef_vect)#np.sum(est_Sxys*this_coef_vect)
                    est_corrs[fcnt]=est_Sxy/np.sqrt(Gzz[0]*psd_mat[this_patch[thislabw],fcnt])
                est_corr=np.mean(est_corrs)
                act_corr=con_mat_mean[cnt1,cnt2]
                if np.abs(act_corr-est_corr)/act_corr<0.1:
                    con_mat_modified[cnt1,cnt2]=0#np.sqrt(psd_mat[cnt2,fcnt]*Gzz[0])#
                else:
                    con_mat_modified[cnt1,cnt2]=act_corr
    con_mat_final=np.minimum(con_mat_modified,np.transpose(con_mat_modified))
#    for cnt1 in range(len(labels)): #seed label
#        for cnt2 in range(len(labels)): #destination labels
#            if cnt1<cnt2:
#               con_mat_final[cnt1,cnt2]=0 
    return con_mat_final



def estimate_actual_signal(labels, Rmat, labels_ts ):
    sRmat=np.divide(Rmat,np.max(Rmat,axis=0))  
    label_patches=np.where(sRmat>0.25)   
    labels_ts_modified=copy.deepcopy(labels_ts) 
    for ecnt in range(len(labels_ts)):
        for cnt2 in range(len(labels)): #destination labels
            this_patch=label_patches[0][label_patches[1]==cnt2]
            thislabw=np.where(this_patch==cnt2)[0]
            coef_mat=Rmat[this_patch,:][:,this_patch].transpose()
            sig_vect=labels_ts_modified[ecnt][this_patch,:]
            est_sig=np.linalg.lstsq(coef_mat,sig_vect)
            labels_ts_modified[ecnt][cnt2,:]=est_sig
            for tcnt in range(labels_ts[ecnt].shape[1]):
                sig_vect=labels_ts_modified[ecnt][this_patch,tcnt]
                est_sig=np.linalg.lstsq(coef_mat,sig_vect)
                labels_ts_modified[ecnt][cnt2,tcnt]=est_sig
    return labels_ts_modified