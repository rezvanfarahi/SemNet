clear
clc
% close all

in_path='/imaging/rf02/TypLexMEG/';
label_path='/imaging/rf02/TypLexMEG/parcellation/labels_and_CTFs/';

list_all = {'meg10_0378/101209/',
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
    };

% load(['/imaging/rf02/TypLexMEG/','3conn_words_',seed,'_',band,'.mat')
% load('/imaging/rf02/TypLexMEG/3conn_nonwords.mat')
% Xc=Xw; Xa=Xnw;
%Xc an Xa: 20484x2x3x17; 3 conn methods: mi, dwpli, coh
load('/imaging/rf02/TypLexMEG/meg10_0378/101209/MIg_14ROIs_2_abstract.mat')

exp_mtx_c=nan(17,2,4);
C_hmuch_c=nan(17,2,4);
C_which_c=nan(17,2,4);
exp_mtx_a=nan(17,2,4);
C_hmuch_a=nan(17,2,4);
C_which_a=nan(17,2,4);
C_c=nan(17,2,2,4);%nsub,ntimes,nconn,nbands
C_a=nan(17,2,2,4);
exp_mtx_bands_c=nan(17,2);
C_hmuch_bands_c=nan(17,2);
C_which_bands_c=nan(17,2);
exp_mtx_bands_a=nan(17,2);
C_hmuch_bands_a=nan(17,2);
C_which_bands_a=nan(17,2);
C_bands_c=nan(17,2,4);
C_bands_a=nan(17,2,4);
ntimes=2;
nbands=4;
nrois=12;
SC_pca_a=zeros(nrois,nrois,2,4);
SC_pca_alla=zeros(nrois,nrois,2,4,17);
SC_pca_bands_a=zeros(nrois,nrois,2);
SC_pca_allab=zeros(nrois,nrois,2,17);
whichpc=1;
whichpcb=1;
whichconn=1;%1 co, 2 ppc, 3 mi
for sub_no=1:length(list_all)
    sub_no
    load(['/imaging/rf02/TypLexMEG/',list_all{sub_no},'ppc_14ROIs_2_abstract.mat'])
    load(['/imaging/rf02/TypLexMEG/',list_all{sub_no},'coh_14ROIs_2_abstract.mat'])
    for time_no=1:ntimes
        for band_no=1:nbands
            
            
            switch band_no
                case 1
                    load(['/imaging/rf02/TypLexMEG/',list_all{sub_no},'MIt_14ROIs_2_abstract.mat'])
                    this_mi=abs(MIt_final(:,:,time_no)); this_mi=this_mi+this_mi';                                        
                case 2
                    load(['/imaging/rf02/TypLexMEG/',list_all{sub_no},'MIa_14ROIs_2_abstract.mat'])
                    this_mi=abs(MIa_final(:,:,time_no)); this_mi=this_mi+this_mi';
                case 3
                    load(['/imaging/rf02/TypLexMEG/',list_all{sub_no},'MIb_14ROIs_2_abstract.mat'])
                    this_mi=abs(MIb_final(:,:,time_no)); this_mi=this_mi+this_mi';
                case 4
                    load(['/imaging/rf02/TypLexMEG/',list_all{sub_no},'MIg_14ROIs_2_abstract.mat'])
                    this_mi=abs(MIg_final(:,:,time_no)); this_mi=this_mi+this_mi';
            end
            this_mi=[this_mi(1:10,:);mean(this_mi([11,13,15],:));mean(this_mi([12,14,16],:))];
            this_mi=[this_mi(:,1:10),mean(this_mi(:,[11,13,15]),2),mean(this_mi(:,[12,14,16]),2)];
            this_mi=this_mi-diag(diag(this_mi));
            this_mi=tril(this_mi,-1);
            this_mi=this_mi(this_mi>0);
            this_ppc=abs(ppc(:,:,band_no,time_no)); 
            this_ppc=[this_ppc(1:10,:);mean(this_ppc([11,13,15],:));mean(this_ppc([12,14,16],:))];
            this_ppc=[this_ppc(:,1:10),mean(this_ppc(:,[11,13,15]),2),mean(this_ppc(:,[12,14,16]),2)];
            this_ppc=this_ppc-diag(diag(this_ppc));
            this_ppc=tril(this_ppc,-1);
            this_ppc=this_ppc(this_ppc>0);
            
            this_coh=abs(coh(:,:,band_no,time_no)); 
            this_coh=[this_coh(1:10,:);mean(this_coh([11,13,15],:));mean(this_coh([12,14,16],:))];
            this_coh=[this_coh(:,1:10),mean(this_coh(:,[11,13,15]),2),mean(this_coh(:,[12,14,16]),2)];
            this_coh=this_coh-diag(diag(this_coh));
            this_coh=tril(this_coh,-1);
            this_coh=this_coh(this_coh>0);
            
            X1=this_mi;%;%[this_coh,this_ppc];%,this_mi];%X1;%-repmat(mean(X1),size(X1,1),1);
            X=X1-repmat(mean(X1),size(X1,1),1);
            [coeff,score,latent,tsquared,explained,mu]=pca(X);
            aa=tril(ones(size(SC_pca_a,1),size(SC_pca_a,2)),-1);
            aa(aa>0)=score(:,whichpc);
            SC_pca_a(:,:,time_no,band_no)=aa;
            
            exp_mtx_a(sub_no,time_no,band_no)=explained(whichpc);
            [C_hmuch_a(sub_no,time_no,band_no),C_which_a(sub_no,time_no,band_no)]=max(corr(X,score(:,whichpc)));
            C_a(sub_no,time_no,:,band_no)=corr(X,score(:,whichpc));
            if band_no==1
                X1_bands=zeros(size(score,1),nbands);
            end
            X1_bands(:,band_no)=score(:,whichpc);%X1(:,whichconn);%
            
        end
        X=X1_bands-repmat(mean(X1_bands),size(X1_bands,1),1);
        [coeff,score,latent,tsquared,explained,mu]=pca(X);
        aa=tril(ones(size(SC_pca_bands_a,1),size(SC_pca_bands_a,2)),-1);
        aa(aa>0)=score(:,whichpcb);
        SC_pca_bands_a(:,:,time_no)=aa;
        exp_mtx_bands_a(sub_no,time_no)=explained(whichpcb);
        [C_hmuch_bands_a(sub_no,time_no),C_which_bands_a(sub_no,time_no)]=max(corr(X,score(:,whichpcb)));
        C_bands_a(sub_no,time_no,:)=corr(X,score(:,whichpcb));
    end
    out_path=['/imaging/rf02/TypLexMEG/',list_all{sub_no},'pca_14ROIs_2_abstract.mat'];
    save(out_path,'SC_pca_a')
    SC_pca_alla(:,:,:,:,sub_no)=SC_pca_a;
    SC_pca_allab(:,:,:,sub_no)=SC_pca_bands_a;
end

SC_pca_c=zeros(nrois,nrois,2,4);
SC_pca_allc=zeros(nrois,nrois,2,4,17);
SC_pca_bands_c=zeros(nrois,nrois,2);
SC_pca_allcb=zeros(nrois,nrois,2,17);
for sub_no=1:length(list_all)
    sub_no
    load(['/imaging/rf02/TypLexMEG/',list_all{sub_no},'ppc_14ROIs_2_concrete.mat'])
    load(['/imaging/rf02/TypLexMEG/',list_all{sub_no},'coh_14ROIs_2_concrete.mat'])
    for time_no=1:ntimes
        for band_no=1:nbands
            
            
            switch band_no
                case 1
                    load(['/imaging/rf02/TypLexMEG/',list_all{sub_no},'MIt_14ROIs_2_concrete.mat'])
                    this_mi=abs(MIt_final(:,:,time_no)); this_mi=this_mi+this_mi';
                case 2
                    load(['/imaging/rf02/TypLexMEG/',list_all{sub_no},'MIa_14ROIs_2_concrete.mat'])
                    this_mi=abs(MIa_final(:,:,time_no)); this_mi=this_mi+this_mi';
                case 3
                    load(['/imaging/rf02/TypLexMEG/',list_all{sub_no},'MIb_14ROIs_2_concrete.mat'])
                    this_mi=abs(MIb_final(:,:,time_no)); this_mi=this_mi+this_mi';
                case 4
                    load(['/imaging/rf02/TypLexMEG/',list_all{sub_no},'MIg_14ROIs_2_concrete.mat'])
                    this_mi=abs(MIg_final(:,:,time_no)); this_mi=this_mi+this_mi';
            end
            this_mi=[this_mi(1:10,:);mean(this_mi([11,13,15],:));mean(this_mi([12,14,16],:))];
            this_mi=[this_mi(:,1:10),mean(this_mi(:,[11,13,15]),2),mean(this_mi(:,[12,14,16]),2)];
            this_mi=this_mi-diag(diag(this_mi));
            this_mi=tril(this_mi,-1);
            this_mi=this_mi(this_mi>0);
            this_ppc=abs(ppc(:,:,band_no,time_no)); 
            this_ppc=[this_ppc(1:10,:);mean(this_ppc([11,13,15],:));mean(this_ppc([12,14,16],:))];
            this_ppc=[this_ppc(:,1:10),mean(this_ppc(:,[11,13,15]),2),mean(this_ppc(:,[12,14,16]),2)];
            this_ppc=this_ppc-diag(diag(this_ppc));
            this_ppc=tril(this_ppc,-1);
            this_ppc=this_ppc(this_ppc>0);            
            this_coh=abs(coh(:,:,band_no,time_no)); 
            this_coh=[this_coh(1:10,:);mean(this_coh([11,13,15],:));mean(this_coh([12,14,16],:))];
            this_coh=[this_coh(:,1:10),mean(this_coh(:,[11,13,15]),2),mean(this_coh(:,[12,14,16]),2)];
            this_coh=this_coh-diag(diag(this_coh));
            this_coh=tril(this_coh,-1);
            this_coh=this_coh(this_coh>0);
            X1=this_mi;%[this_coh,this_ppc];%,this_mi];%X1;%-repmat(mean(X1),size(X1,1),1);
            X=X1-repmat(mean(X1),size(X1,1),1);
            [coeff,score,latent,tsquared,explained,mu]=pca(X);
            aa=tril(ones(size(SC_pca_c,1),size(SC_pca_c,2)),-1);
            aa(aa>0)=score(:,whichpc);
            SC_pca_c(:,:,time_no,band_no)=aa;
            exp_mtx_c(sub_no,time_no,band_no)=explained(whichpc);
            [C_hmuch_c(sub_no,time_no,band_no),C_which_c(sub_no,time_no,band_no)]=max(corr(X,score(:,whichpc)));
            C_c(sub_no,time_no,:,band_no)=corr(X,score(:,whichpc));
            if band_no==1
                X1_bands=zeros(size(score,1),nbands);
            end
            X1_bands(:,band_no)=score(:,whichpc);%X1(:,whichconn);%
            
            
        end
        X=X1_bands-repmat(mean(X1_bands),size(X1_bands,1),1);
        [coeff,score,latent,tsquared,explained,mu]=pca(X);
        aa=tril(ones(size(SC_pca_bands_c,1),size(SC_pca_bands_c,2)),-1);
        aa(aa>0)=score(:,whichpcb);
        SC_pca_bands_c(:,:,time_no)=aa;
        exp_mtx_bands_c(sub_no,time_no)=explained(whichpcb);
        [C_hmuch_bands_c(sub_no,time_no),C_which_bands_c(sub_no,time_no)]=max(corr(X,score(:,whichpcb)));
        C_bands_c(sub_no,time_no,:)=corr(X,score(:,whichpcb));
    end
    out_path=['/imaging/rf02/TypLexMEG/',list_all{sub_no},'pca_14ROIs_2_concrete.mat'];
    save(out_path,'SC_pca_c')
    SC_pca_allc(:,:,:,:,sub_no)=SC_pca_c;
    SC_pca_allcb(:,:,:,sub_no)=SC_pca_bands_c;
end

% coh_perc=(sum(C_which_c==3)+sum(C_which_a==3))./(sum(C_which_c>0)+sum(C_which_a>0));
% coh_perc=mean(coh_perc(:));


% clear
% clc
% seeds={'lhATL','lhpIFG','lhMTL','lhSMG'};
% bands={'Theta','Alpha','Beta','Gamma'};
% % load(['/imaging/rf02/TypLexMEG/','3conn_words_',seed,'_',band,'.mat')
% % load('/imaging/rf02/TypLexMEG/3conn_nonwords.mat')
% % Xc=Xw; Xa=Xnw;
% %Xc an Xa: 20484x2x3x17; 3 conn methods: mi, dwpli, coh
% exp_mtx_c=nan(20484,2,4,4);
% C_hmuch_c=nan(20484,2,4,4);
% C_which_c=nan(20484,2,4,4);
% exp_mtx_a=nan(20484,2,4,4);
% C_hmuch_a=nan(20484,2,4,4);
% C_which_a=nan(20484,2,4,4);
% C_c=nan(20484,2,3,4,4);
% C_a=nan(20484,2,3,4,4);
% for seed_no=1:length(seeds)
%     seed=seeds{seed_no};
%     for band_no=1:length(bands)
%         band=bands{band_no};
%         load(['/imaging/rf02/TypLexMEG/','3conn_cncrt_',seed,'_',band,'.mat'])
%         load(['/imaging/rf02/TypLexMEG/','3conn_abs_',seed,'_',band,'.mat'])
%         seed_no
%         band_no
%         for vert_no=1:size(Xc,1)
%             for time_no=1:size(Xc,2)
%                 %concrete
%                 X=squeeze(Xc(vert_no,time_no,:,:)-Xa(vert_no,time_no,:,:))';
%                 [coeff,score,latent,tsquared,explained,mu]=pca(X);
%                 exp_mtx_c(vert_no,time_no,band_no,seed_no)=explained(1);
%                 [C_hmuch_c(vert_no,time_no,band_no,seed_no),C_which_c(vert_no,time_no,band_no,seed_no)]=max(corr(X,score(:,1)));
%                 C_c(vert_no,time_no,:,band_no,seed_no)=corr(X,score(:,1));
%                 %abstract
% %                 X=squeeze(Xa(vert_no,time_no,:,:))';
% %                 [coeff,score,latent,tsquared,explained,mu]=pca(X);
% %                 exp_mtx_a(vert_no,time_no,band_no,seed_no)=explained(1);
% %                 [C_hmuch_a(vert_no,time_no,band_no,seed_no),C_which_a(vert_no,time_no,band_no,seed_no)]=max(corr(X,score(:,1)));
% %                 C_a(vert_no,time_no,:,band_no,seed_no)=corr(X,score(:,1));
%             end
%         end
%     end
% end
% E=mean(exp_mtx_a,4);
E=squeeze(mean(exp_mtx_a,1));
E=squeeze(mean(E,1));

% E1=mean(C_a,4);
E1=squeeze(mean(C_a,1));
E1=squeeze(mean(E1,1));
PCA_SLa=[E;E1];
E=squeeze(mean(exp_mtx_c,1));
E=squeeze(mean(E,1));

% E1=mean(C_a,4);
E1=squeeze(mean(C_c,1));
E1=squeeze(mean(E1,1));
PCA_SLc=[E;E1];
cd('/home/rf02/rezvan/test1/figures')
% save('PCA_acrosssubjs_SL.mat','PCA_SL')
% coh_perc=(sum(C_which_c==2)+sum(C_which_a==2))./(sum(C_which_c>0)+sum(C_which_a>0));
% coh_perc=mean(coh_perc(:));
%%% PCA bands
Eb=squeeze(mean(exp_mtx_bands_a,1));
Eb=squeeze(mean(Eb));

% E1=mean(C_a,4);
E1b=squeeze(mean(C_bands_a,1));
E1b=squeeze(mean(E1b,1));
PCA_SLab=[Eb,E1b];
Eb=squeeze(mean(exp_mtx_bands_c,1));
Eb=squeeze(mean(Eb));

% E1=mean(C_a,4);
E1b=squeeze(mean(C_bands_c,1));
E1b=squeeze(mean(E1b,1));
PCA_SLcb=[Eb,E1b];

SC_pca_alls=SC_pca_allcb-SC_pca_allab;
pb=zeros(nrois,nrois,2);
tb=zeros(nrois,nrois,2);
for c1=1:size(SC_pca_alls,1)
    for c2=1:size(SC_pca_alls,2)
        for c3=1:size(SC_pca_alls,3)
            
                this_vect=SC_pca_alls(c1,c2,c3,:);
                [h,pb(c1,c2,c3),ci,stats]=ttest(this_vect);
                tb(c1,c2,c3)=stats.tstat;
            
        end
    end
end
figure,
cnt=0;
for jj=1:2
cnt=cnt+1;
subplot(1,2,cnt); imagesc(squeeze(tb(:,:,jj)))
end

colormap('cool')
colorbar()

SC_pca_alls=SC_pca_allc-SC_pca_alla;
p=zeros(nrois,nrois,2,4);
t=zeros(nrois,nrois,2,4);
for c1=1:size(SC_pca_alls,1)
    for c2=1:size(SC_pca_alls,2)
        for c3=1:size(SC_pca_alls,3)
            for c4=1:size(SC_pca_alls,4)
                this_vect=SC_pca_alls(c1,c2,c3,c4,:);
                [h,p(c1,c2,c3,c4),ci,stats]=ttest(this_vect);
                t(c1,c2,c3,c4)=stats.tstat;
            end
        end
    end
end
figure,
cnt=0;
for ii=1:4
for jj=1:2
cnt=cnt+1;
subplot(2,4,cnt); imagesc(squeeze(t(:,:,jj,ii)))
end
end
colormap('cool')
colorbar()

% out_path='/imaging/rf02/TypLexMEG/icaanalysis_results/ROI/concrete_coh_ppc_1stpc_bands_1stpc.mat';
% save(out_path,'SC_pca_allcb')
% out_path='/imaging/rf02/TypLexMEG/icaanalysis_results/ROI/abstract_coh_ppc_1stpc_bands_1stpc.mat';
% save(out_path,'SC_pca_allab')
% 
% out_path='/imaging/rf02/TypLexMEG/icaanalysis_results/ROI/concrete_coh_ppc_1stpc_bands_all.mat';
% save(out_path,'SC_pca_allc')
% out_path='/imaging/rf02/TypLexMEG/icaanalysis_results/ROI/abstract_coh_ppc_1stpc_bands_all.mat';
% save(out_path,'SC_pca_alla')

out_path='/imaging/rf02/TypLexMEG/icaanalysis_results/ROI/concrete_mi_bands_1stpc.mat';
save(out_path,'SC_pca_allcb')
out_path='/imaging/rf02/TypLexMEG/icaanalysis_results/ROI/abstract_mi_bands_1stpc.mat';
save(out_path,'SC_pca_allab')

out_path='/imaging/rf02/TypLexMEG/icaanalysis_results/ROI/concrete_mi_bands_all.mat';
save(out_path,'SC_pca_allc')
out_path='/imaging/rf02/TypLexMEG/icaanalysis_results/ROI/abstract_mi_bands_all.mat';
save(out_path,'SC_pca_alla')
