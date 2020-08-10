clear
clc
bands={'theta','alpha','beta','gamma'};
% load(['/imaging/rf02/TypLexMEG/','3conn_words_',seed,'_',band,'.mat')
% load('/imaging/rf02/TypLexMEG/3conn_nonwords.mat')
% Xc=Xw; Xa=Xnw;
%Xc an Xa: 20484x2x3x17; 3 conn methods: mi, dwpli, coh#

%% over verts- bands
ntimes=5;
nsubs=17;
nverts=20484;
% which_pc=4;
nbands=length(bands);
npcs=nbands;
exp_mtx_c=nan(nsubs,ntimes,npcs);
C_hmuch_c=nan(nsubs,ntimes,npcs);
C_which_c=nan(nsubs,ntimes,npcs);
exp_mtx_a=nan(nsubs,ntimes,npcs);
C_hmuch_a=nan(nsubs,ntimes,npcs);
C_which_a=nan(nsubs,ntimes,npcs);
exp_mtx_s=nan(nsubs,ntimes,npcs);
C_hmuch_s=nan(nsubs,ntimes,npcs);
C_which_s=nan(nsubs,ntimes,npcs);
exp_mtx_b=nan(nsubs,ntimes,npcs);
C_hmuch_b=nan(nsubs,ntimes,npcs);
C_which_b=nan(nsubs,ntimes,npcs);
C_c=nan(nsubs,ntimes,nbands,npcs);
C_a=nan(nsubs,ntimes,nbands,npcs);
C_s=nan(nsubs,ntimes,nbands,npcs);
C_b=nan(nsubs,ntimes,nbands,npcs);

Xc_score=zeros(nverts,ntimes,nsubs,npcs);
Xa_score=zeros(nverts,ntimes,nsubs,npcs);
Xs_score=zeros(nverts,ntimes,nsubs,npcs);
Xb_score=zeros(2*nverts,ntimes,nsubs,npcs);

Xc_grand=zeros(nverts,ntimes,nsubs,nbands);
Xa_grand=zeros(nverts,ntimes,nsubs,nbands);
Xs_grand=zeros(nverts,ntimes,nsubs,nbands);
Xb_grand=zeros(2*nverts,ntimes,nsubs,nbands);
for band_no=1:nbands
    band=bands{band_no};
    load(['/imaging/rf02/TypLexMEG/','TF_orig_cncrt_',band,'.mat'])
    load(['/imaging/rf02/TypLexMEG/','TF_orig_abs_',band,'.mat'])
    
    Xc=abs(Xc-1);
    Xa=abs(Xa-1);
    Xc_grand(:,:,:,band_no)=Xc;
    Xa_grand(:,:,:,band_no)=Xa;
    Xs_grand(:,:,:,band_no)=Xc-Xa;
    Xb_grand(:,:,:,band_no)=cat(1,Xc,Xa);
end
for sub_no=1:nsubs
    sub_no
    for time_no=1:ntimes
        for which_pc=1:nbands
        
        %concrete
        X1=squeeze(Xc_grand(:,time_no,sub_no,:));%-Xa(:,time_no,:,sub_no));
        X=X1;%-repmat(mean(X1),size(X1,1),1);
        [coeff,score,latent,tsquared,explained,mu]=pca(X,'Centered',true);%'VariableWeights','variance','VariableWeights','variance',
        exp_mtx_c(sub_no,time_no,:)=explained;%(which_pc);
        [C_hmuch_c(sub_no,time_no,which_pc),C_which_c(sub_no,time_no,which_pc)]=max(corr(X,score(:,which_pc)));
        Xc_score(:,time_no,sub_no,:)=score;%-repmat(min(score),size(score,1),1);
        C_c(sub_no,time_no,:,:)= corr(X,score);%coeff;%corr(X,score(:,which_pc));inv(diag(std(X)))* coeff;
        %abstract
        X1=squeeze(Xa_grand(:,time_no,sub_no,:));
        X=X1;%-repmat(mean(X1),size(X1,1),1);
        [coeff,score,latent,tsquared,explained,mu]=pca(X,'Centered',true);
        exp_mtx_a(sub_no,time_no,:)=explained;%(which_pc);
        [C_hmuch_a(sub_no,time_no,which_pc),C_which_a(sub_no,time_no,which_pc)]=max(corr(X,score(:,which_pc)));
        Xa_score(:,time_no,sub_no,:)=score;%-repmat(min(score),size(score,1),1);
        C_a(sub_no,time_no,:,:)= corr(X,score);%coeff;%corr(X,score(:,which_pc));
        %subtract
        X1=squeeze(Xs_grand(:,time_no,sub_no,:));
        X=X1;%-repmat(mean(X1),size(X1,1),1);
        [coeff,score,latent,tsquared,explained,mu]=pca(X,'Centered',true);
        exp_mtx_s(sub_no,time_no,:)=explained;%(which_pc);
        [C_hmuch_s(sub_no,time_no,which_pc),C_which_s(sub_no,time_no,which_pc)]=max(corr(X,score(:,which_pc)));
        Xs_score(:,time_no,sub_no,:)=score;%-repmat(min(score),size(score,1),1);
        C_s(sub_no,time_no,:,:)= corr(X,score);%coeff;%corr(X,score(:,which_pc));
        
        %pool both
        X1=squeeze(Xb_grand(:,time_no,sub_no,:));
        X=X1;%-repmat(mean(X1),size(X1,1),1);
        [coeff,score,latent,tsquared,explained,mu]=pca(X,'Centered',true);
        exp_mtx_b(sub_no,time_no,:)=explained;%(which_pc);
        [C_hmuch_s(sub_no,time_no,which_pc),C_which_s(sub_no,time_no,which_pc)]=max(corr(X,score(:,which_pc)));
        Xb_score(:,time_no,sub_no,:)=score;%-repmat(min(score),size(score,1),1);
        C_b(sub_no,time_no,:,:)=corr(X,score);%coeff;%corr(X,score(:,which_pc));
        end 
    end
end


E=squeeze(mean(exp_mtx_a,1));
E=squeeze(mean(E,1));

E1=squeeze(mean(C_a,1));
E1=squeeze(mean(E1,1));
PCA_SLa=[E;E1];

E=squeeze(mean(exp_mtx_c,1));
E=squeeze(mean(E,1));

E1=squeeze(mean(C_c,1));
E1=squeeze(mean(E1,1));
PCA_SLc=[E;E1];

E=mean(exp_mtx_s,4);
E=squeeze(mean(E,1));
E=squeeze(mean(E,1));

E1=squeeze(mean(C_s,1));
E1=squeeze(mean(E1,1));
PCA_SLs=[E;E1];

E=mean(exp_mtx_b,4);
E=squeeze(mean(E,1));
E=squeeze(mean(E,1));

E1=squeeze(mean(C_b,1));
E1=squeeze(mean(E1,1));
PCA_SLb=[E;E1];

PCA_SL_ca=[PCA_SLc;PCA_SLa;PCA_SLs];
PCA_SL_sb=[PCA_SLb;PCA_SLs];
% out_path=['/imaging/rf02/TypLexMEG/icaanalysis_results/tf_pca/WB_PCA_acrossverts_TL_allbands.mat'];
% save(out_path,'PCA_SL_ca')
% out_path=['/imaging/rf02/TypLexMEG/icaanalysis_results/tf_pca/WB_PCA_acrossverts_TL_allbands_Xcscore.mat'];
% save(out_path,'Xc_score')
% 
% out_path=['/imaging/rf02/TypLexMEG/icaanalysis_results/tf_pca/WB_PCA_acrossverts_TL_allbands_Xascore.mat'];
% save(out_path,'Xa_score')
% 
% out_path=['/imaging/rf02/TypLexMEG/icaanalysis_results/tf_pca/WB_PCA_acrossverts_TL_allbands_Xsscore.mat'];
% save(out_path,'Xs_score')
% 
% out_path=['/imaging/rf02/TypLexMEG/icaanalysis_results/tf_pca/WB_PCA_acrossverts_TL_allbands_Xbscore.mat'];
% save(out_path,'Xb_score')

out_path=['/imaging/rf02/TypLexMEG/icaanalysis_results/tf_pca/WB_PCA_acrossverts_TL_allbands_Xcgrand.mat'];
save(out_path,'Xc_grand')

out_path=['/imaging/rf02/TypLexMEG/icaanalysis_results/tf_pca/WB_PCA_acrossverts_TL_allbands_Xagrand.mat'];
save(out_path,'Xa_grand')

out_path=['/imaging/rf02/TypLexMEG/icaanalysis_results/tf_pca/WB_PCA_acrossverts_TL_allbands_Xsgrand.mat'];
save(out_path,'Xs_grand')

out_path=['/imaging/rf02/TypLexMEG/icaanalysis_results/tf_pca/WB_PCA_acrossverts_TL_allbands_Xbgrand.mat'];
save(out_path,'Xb_grand')



% xx1=squeeze(mean(Xc_grand,2));%over time
% xx1=squeeze(mean(xx1,2));%over subs
% xx2=squeeze(mean(Xa_grand,2));%over time
% xx2=squeeze(mean(xx2,2));%over subs
% xxborig=corr(xx1,xx2);
%  subplot(2,1,1);imagesc(xxborig)%figure,
% colormap('hot')
% colorbar()
Xborig=zeros(nbands,nbands,nsubs);
for jj=1:size(Xc_grand,3) %subj
    xxx=[];
for ii = 1:size(Xc_grand,2) %time
   
        
xx1=squeeze(Xc_grand(:,ii,jj,:));
xx2=squeeze(Xa_grand(:,ii,jj,:));%over time
xxborig=corr(xx1,xx2);
xxx=cat(3,xxx,xxborig);
end
Xborig(:,:,jj)=mean(xxx,3);
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
        
Xbpca=zeros(nbands,nbands,nsubs);
xcs=Xb_score(1:20484,:,:,:);
xas=Xb_score(20485:end,:,:,:);

for jj=1:size(Xc_grand,3) %subj
    xxx=[];
for ii = 1:size(Xc_grand,2) %time
   
        
xx1=squeeze(xcs(:,ii,jj,:));
xx2=squeeze(xas(:,ii,jj,:));%over time
xxborig=corr(xx1,xx2);
xxx=cat(3,xxx,xxborig);
end
Xbpca(:,:,jj)=mean(xxx,3);
end
 
Xbpca_pval=zeros(size(Xbpca,1),size(Xbpca,2));
for p1cnt=1:size(Xbpca_pval,1)
    for p2cnt = 1:size(Xbpca_pval,2)
        [h,Xbpca_pval(p1cnt,p2cnt),ci,stats]= ttest(squeeze(Xbpca(p1cnt,p2cnt,:)));
    end
end
corrmatMpca=mean(Xbpca,3);


% xcs=Xb_score(1:20484,:,:,:);
% xx1=squeeze(mean(xcs,2));%over time
% xx1=squeeze(mean(xx1,2));%over subs
% xas=Xb_score(20485:end,:,:,:);
% xx2=squeeze(mean(xas,2));%over time
% xx2=squeeze(mean(xx2,2));%over subs
% xxbpc=corr(xx1,xx2);
subplot(2,1,2); imagesc(corrmatMpca)
colormap('hot')
colorbar()