clear
clc
seeds={'lhATL','lhpIFG','lhMTL','lhAG'};%'rhATL',
bands={'Theta','Alpha','Beta','Gamma'};
connms={'mi','ppc','coh'};
% load(['/imaging/rf02/TypLexMEG/','3conn_words_',seed,'_',band,'.mat')
% load('/imaging/rf02/TypLexMEG/3conn_nonwords.mat')
% Xc=Xw; Xa=Xnw;
%Xc an Xa: 20484x2x3x17; 3 conn methods: mi, dwpli, coh#
%% over verts- bands
ntimes=2;
nsubs=17;
nverts=20484;
% which_pc=4;
nbands=length(bands);
npcs=nbands;
nseeds=length(seeds);


Xc_score=zeros(nverts,ntimes,nsubs,npcs);
Xa_score=zeros(nverts,ntimes,nsubs,npcs);
Xs_score=zeros(nverts,ntimes,nsubs,npcs);
Xb_score=zeros(2*nverts,ntimes,nsubs,npcs);

Xc_grand=zeros(nverts,ntimes,nsubs,nbands,nseeds);
Xa_grand=zeros(nverts,ntimes,nsubs,nbands,nseeds);
Xs_grand=zeros(nverts,ntimes,nsubs,nbands,nseeds);
Xb_grand=zeros(2*nverts,ntimes,nsubs,nbands,nseeds);

Xc_gscore=zeros(nverts,ntimes,nsubs,nbands,nseeds);
Xa_gscore=zeros(nverts,ntimes,nsubs,nbands,nseeds);
Xs_gscore=zeros(nverts,ntimes,nsubs,nbands,nseeds);
Xb_gscore=zeros(2*nverts,ntimes,nsubs,nbands,nseeds);

which_conn=1;
for seed_no=1:nseeds
    seed=seeds{seed_no};
    for band_no=1:nbands
        band=bands{band_no};
        load(['/imaging/rf02/TypLexMEG/','3conn_newag_orig_cncrt_',seed,'_',band,'.mat'])
        load(['/imaging/rf02/TypLexMEG/','3conn_newag_orig_abs_',seed,'_',band,'.mat'])
        if which_conn==1
            for xc1=1:size(Xc,2)
                for xc2=1:size(Xc,4)
                    [cc,ccw]=sort(Xc(:,xc1,which_conn,xc2),'descend');
                    [aa,aaw]=sort(Xa(:,xc1,which_conn,xc2),'descend');
                    Xc(ccw(1),xc1,which_conn,xc2)=cc(2);
                    Xa(aaw(1),xc1,which_conn,xc2)=aa(2);
                    
                end
            end
        end
        
        Xc_grand(:,:,:,band_no,seed_no)=squeeze(Xc(:,:,which_conn,:));
        Xa_grand(:,:,:,band_no,seed_no)=squeeze(Xa(:,:,which_conn,:));
        
        Xs_grand(:,:,:,band_no,seed_no)=squeeze(Xc(:,:,which_conn,:))-squeeze(Xa(:,:,which_conn,:));
        Xb_grand(:,:,:,band_no,seed_no)=cat(1,squeeze(Xc(:,:,which_conn,:)),squeeze(Xa(:,:,which_conn,:)));
    end
    for sub_no=1:nsubs
        sub_no
        for time_no=1:ntimes
            for which_pc=1:nbands
                
                %concrete
                X1=squeeze(Xc_grand(:,time_no,sub_no,:,seed_no));%-Xa(:,time_no,:,sub_no));
                X=X1;%-repmat(mean(X1),size(X1,1),1);
                [coeff,scorec,latent,tsquared,explained,mu]=pca(X,'Centered',true);%,'Centered',true);
                Xc_score(:,time_no,sub_no,:)=scorec;%-repmat(min(scorec),size(scorec,1),1);
                Xc_gscore(:,time_no,sub_no,:,seed_no)=scorec;
                %abstract
                X1=squeeze(Xa_grand(:,time_no,sub_no,:,seed_no));
                X=X1;%-repmat(mean(X1),size(X1,1),1);
                [coeff,scorea,latent,tsquared,explained,mu]=pca(X,'Centered',true);%,'Centered',true);
                Xa_score(:,time_no,sub_no,:)=scorea;%-repmat(min(scorea),size(scorea,1),1);
                Xa_gscore(:,time_no,sub_no,:,seed_no)=scorea;
                %subtract
                X1=squeeze(Xs_grand(:,time_no,sub_no,:,seed_no));
                X=X1;%-repmat(mean(X1),size(X1,1),1);
                [coeff,score,latent,tsquared,explained,mu]=pca(X,'Centered',true);%,'Centered',true);
                Xs_score(:,time_no,sub_no,:)=score;%-repmat(min(score),size(score,1),1);
                Xs_gscore(:,time_no,sub_no,:,seed_no)=score;
                %both pooled
                X1=squeeze(Xb_grand(:,time_no,sub_no,:,seed_no));
                X=X1;%-repmat(mean(X1),size(X1,1),1);
                [coeff,score,latent,tsquared,explained,mu]=pca(X,'Centered',true);%,'Centered',true);
                Xb_score(:,time_no,sub_no,:)=score;%-repmat(min(score),size(score,1),1);
                Xb_gscore(:,time_no,sub_no,:,seed_no)=score;
            end
        end
    end
    
%     out_path=['/imaging/rf02/TypLexMEG/icaanalysis_results/conn_pca/WB_ncPCA_',connms{which_conn}, '_acrossverts_TL_allbands_Xcscore_',seed,'.mat'];
%     save(out_path,'Xc_score')
%     
%     out_path=['/imaging/rf02/TypLexMEG/icaanalysis_results/conn_pca/WB_ncPCA_',connms{which_conn}, '_acrossverts_TL_allbands_Xascore_',seed,'.mat'];
%     save(out_path,'Xa_score')
%     
%     out_path=['/imaging/rf02/TypLexMEG/icaanalysis_results/conn_pca/WB_ncPCA_',connms{which_conn}, '_acrossverts_TL_allbands_Xsscore_',seed,'.mat'];
%     save(out_path,'Xs_score')
%     out_path=['/imaging/rf02/TypLexMEG/icaanalysis_results/conn_pca/WB_ncPCA_',connms{which_conn}, '_acrossverts_TL_allbands_Xbscore_',seed,'.mat'];
%     save(out_path,'Xb_score')
end
% U, s, V = linalg.svd(stc.data[vertidx, :],
%                                          full_matrices=False)
%                     # determine sign-flip
%                     sign = np.sign(np.dot(U[:, 0], flip))
% 
%                     # use average power in label for scaling
%                     scale = linalg.norm(s) / np.sqrt(len(vertidx))
% 
%                     label_tc[i] = sign * scale * V[0]
Xborig=zeros(nbands,nbands,nsubs);
for jj=1:size(Xc_grand,3) %subj
    xxs=[];
    for kk=1:size(Xc_grand,5)
        xxt=[];
        for ii = 1:size(Xc_grand,2) %time
            xx1=squeeze(Xc_grand(:,ii,jj,:,kk));
            xx2=squeeze(Xa_grand(:,ii,jj,:,kk));%over time
            xxborig=corr(xx1,xx2);
            xxt=cat(3,xxt,xxborig);
        end
        xxs=cat(3,xxs,mean(xxt,3));
    end
    Xborig(:,:,jj)=mean(xxt,3);
end
subplot(2,1,1);imagesc(mean(Xborig,3))%figure,
colormap('hot')
colorbar()
Xborig_pval=zeros(size(Xborig,1),size(Xborig,2));
for p1cnt=1:size(Xborig_pval,1)
    for p2cnt = 1:size(Xborig_pval,2)
        [h,Xborig_pval(p1cnt,p2cnt),ci,stats]= ttest(squeeze(Xborig(p1cnt,p2cnt,:)));
    end
end
corrmatM=mean(Xborig,3);


xcs=Xb_gscore(1:20484,:,:,:,:);
xas=Xb_gscore(20485:end,:,:,:,:);
Xbpca=zeros(nbands,nbands,nsubs);
for jj=1:size(Xb_gscore,3) %subj
    xxs=[];
    for kk=1:size(Xb_gscore,5)
        xxt=[];
        for ii = 1:size(Xb_gscore,2) %time
            xx1=squeeze(xcs(:,ii,jj,:,kk));
            xx2=squeeze(xas(:,ii,jj,:,kk));%over time
            xxbpca=corr(xx1,xx2);
            xxt=cat(3,xxt,xxbpca);
        end
        xxs=cat(3,xxs,mean(xxt,3));
    end
    Xbpca(:,:,jj)=mean(xxt,3);
end
subplot(2,1,2);imagesc(mean(Xbpca,3))%figure,
colormap('hot')
colorbar()
Xbpca_pval=zeros(size(Xbpca,1),size(Xbpca,2));
for p1cnt=1:size(Xbpca_pval,1)
    for p2cnt = 1:size(Xbpca_pval,2)
        [h,Xbpca_pval(p1cnt,p2cnt),ci,stats]= ttest(squeeze(Xbpca(p1cnt,p2cnt,:)));
    end
end
corrmatM=mean(Xbpca,3);



% xx1=squeeze(mean(Xc_gscore,2));%over time
% xx1=squeeze(mean(xx1,4));%over seeds
% xx1=squeeze(mean(xx1,2));%over subs
% xx2=squeeze(mean(Xa_gscore,2));%over time
% xx2=squeeze(mean(xx2,4));%over seeds
% xx2=squeeze(mean(xx2,2));%over subs
% xxb=corr(xx1,xx2);
% figure, imagesc(xxb)
% colormap('hot')
% colorbar()
% 
% xx1=squeeze(mean(Xc_grand,2));%over time
% xx1=squeeze(mean(xx1,4));%over seeds
% xx1=squeeze(mean(xx1,2));%over subs
% xx2=squeeze(mean(Xa_grand,2));%over time
% xx2=squeeze(mean(xx2,4));%over seeds
% xx2=squeeze(mean(xx2,2));%over subs
% xxb=corr(xx1,xx2);
% figure, imagesc(xxb)
% colormap('hot')
% colorbar()