clear
clc
close all
seeds={'lhATL','lhpIFG','lhMTL','lhAG'};%'rhATL',
bands={'Theta','Alpha','Beta','Gamma'};
% load(['/imaging/rf02/TypLexMEG/','3conn_words_',seed,'_',band,'.mat')
% load('/imaging/rf02/TypLexMEG/3conn_nonwords.mat')
% Xc=Xw; Xa=Xnw;
%Xc an Xa: 20484x2x3x17; 3 conn methods: mi, dwpli, coh#
%% over verts - con methods
% connms={'mi','ppc','coh'};%'mi',
% exp_mtx_c=nan(17,2,length(bands),length(seeds));
% C_hmuch_c=nan(17,2,length(bands),length(seeds));
% C_which_c=nan(17,2,length(bands),length(seeds));
% exp_mtx_a=nan(17,2,length(bands),length(seeds));
% C_hmuch_a=nan(17,2,length(bands),length(seeds));
% C_which_a=nan(17,2,length(bands),length(seeds));
% exp_mtx_s=nan(17,2,length(bands),length(seeds));
% C_hmuch_s=nan(17,2,length(bands),length(seeds));
% C_which_s=nan(17,2,length(bands),length(seeds));
% exp_mtx_b=nan(17,2,length(bands),length(seeds));
% C_hmuch_b=nan(17,2,length(bands),length(seeds));
% C_which_b=nan(17,2,length(bands),length(seeds));
% C_c=nan(17,2,length(connms),length(bands),length(seeds));
% C_a=nan(17,2,length(connms),length(bands),length(seeds));
% C_s=nan(17,2,length(connms),length(bands),length(seeds));
% C_b=nan(17,2,length(connms),length(bands),length(seeds));
% which_pc=3;
% for seed_no=1:length(seeds)
%     seed=seeds{seed_no};
%     seed
%     for band_no=1:length(bands)
%         band=bands{band_no};
%         band
%         load(['/imaging/rf02/TypLexMEG/','3conn_newag_orig_cncrt_',seed,'_',band,'.mat'])
%         load(['/imaging/rf02/TypLexMEG/','3conn_newag_orig_abs_',seed,'_',band,'.mat'])
%         for xc1=1:size(Xc,2)
%                 for xc2=1:size(Xc,4)
%                     [cc,ccw]=sort(Xc(:,xc1,1,xc2),'descend');
%                     [aa,aaw]=sort(Xa(:,xc1,1,xc2),'descend');
%                     Xc(ccw(1),xc1,1,xc2)=cc(2);
%                     Xa(aaw(1),xc1,1,xc2)=aa(2);
%                     
%                 end
%         end
%         for xcnt1=1:size(Xc,2)
%             for xcnt2=1:size(Xc,3)
%                 for xcnt3=1:size(Xc,4)
%                     Xc(:,xcnt1,xcnt2,xcnt3)=boxcox1(Xc(:,xcnt1,xcnt2,xcnt3));
%                     Xa(:,xcnt1,xcnt2,xcnt3)=boxcox1(Xa(:,xcnt1,xcnt2,xcnt3));
%                 end
%             end
%         end
%         Xs=Xc-Xa;
%         Xb=cat(1,Xc,Xa);
% 
%         for sub_no=1:size(Xc,4)
% %             sub_no
%             for time_no=1:size(Xc,2)
%                 %concrete
%                 X1=squeeze(Xc(:,time_no,1:end,sub_no));%-Xa(:,time_no,:,sub_no));
%                 X=X1;%-repmat(mean(X1),size(X1,1),1);
%                 [coeff,score,latent,tsquared,explained,mu]=pca(X,'VariableWeights','variance','Centered',true);
%                 exp_mtx_c(sub_no,time_no,band_no,seed_no)=explained(which_pc);
% %                 [C_hmuch_c(sub_no,time_no,band_no,seed_no),C_which_c(sub_no,time_no,band_no,seed_no)]=max(corr(X,score(:,which_pc)));
%                 C_c(sub_no,time_no,:,band_no,seed_no)=corr(X,score(:,which_pc));%inv(diag(std(X)))* coeff(:,which_pc);%
%                 %abstract
%                 X1=squeeze(Xa(:,time_no,1:end,sub_no));
%                 X=X1;%-repmat(mean(X1),size(X1,1),1);
%                 [coeff,score,latent,tsquared,explained,mu]=pca(X,'VariableWeights','variance','Centered',true);
%                 exp_mtx_a(sub_no,time_no,band_no,seed_no)=explained(which_pc);
% %                 [C_hmuch_a(sub_no,time_no,band_no,seed_no),C_which_a(sub_no,time_no,band_no,seed_no)]=max(corr(X,score(:,which_pc)));
%                 C_a(sub_no,time_no,:,band_no,seed_no)=corr(X,score(:,which_pc));%inv(diag(std(X)))* coeff(:,which_pc);%corr(X,score(:,which_pc));
%                 %subtract
%                 X1=squeeze(Xs(:,time_no,1:end,sub_no));
%                 X=X1;%-repmat(mean(X1),size(X1,1),1);
%                 [coeff,score,latent,tsquared,explained,mu]=pca(X,'VariableWeights','variance','Centered',true);
%                 exp_mtx_s(sub_no,time_no,band_no,seed_no)=explained(which_pc);
% %                 [C_hmuch_s(sub_no,time_no,band_no,seed_no),C_which_s(sub_no,time_no,band_no,seed_no)]=max(corr(X,score(:,which_pc)));
%                 C_s(sub_no,time_no,:,band_no,seed_no)=corr(X,score(:,which_pc));%inv(diag(std(X)))* coeff(:,which_pc);%corr(X,score(:,which_pc));
%                 %both
%                 X1=squeeze(Xb(:,time_no,1:end,sub_no));
%                 X=X1;%-repmat(mean(X1),size(X1,1),1);
%                 [coeff,score,latent,tsquared,explained,mu]=pca(X,'VariableWeights','variance','Centered',true);
%                 exp_mtx_b(sub_no,time_no,band_no,seed_no)=explained(which_pc);
% %                 [C_hmuch_s(sub_no,time_no,band_no,seed_no),C_which_s(sub_no,time_no,band_no,seed_no)]=max(corr(X,score(:,which_pc)));
%                 C_b(sub_no,time_no,:,band_no,seed_no)=corr(X,score(:,which_pc));%inv(diag(std(X)))* coeff(:,which_pc);%
%             end
%         end
%     end
% end
% 
% E=mean(exp_mtx_a,4);
% E=squeeze(mean(E,1));
% E=squeeze(mean(E,1));
% E1=mean(C_a,5);
% E1=squeeze(mean(E1,1));
% E1=squeeze(mean(E1,1));
% PCA_SLa=[E;E1];
% 
% E=mean(exp_mtx_c,4);
% E=squeeze(mean(E,1));
% E=squeeze(mean(E,1));
% E1=mean(C_c,5);
% E1=squeeze(mean(E1,1));
% E1=squeeze(mean(E1,1));
% PCA_SLc=[E;E1];
% 
% E=mean(exp_mtx_s,4);
% E=squeeze(mean(E,1));
% E=squeeze(mean(E,1));
% E1=mean(C_s,5);
% E1=squeeze(mean(E1,1));
% E1=squeeze(mean(E1,1));
% PCA_SLs=[E;E1];
% 
% E=mean(exp_mtx_b,4);
% E=squeeze(mean(E,1));
% E=squeeze(mean(E,1));
% E1=mean(C_b,5);
% E1=squeeze(mean(E1,1));
% E1=squeeze(mean(E1,1));
% PCA_SLb=[E;E1];
% PCA_SL_sb=[PCA_SLb;PCA_SLs];
% PCA_SL_ca=[PCA_SLc;PCA_SLa;PCA_SLs];
% out_path=['/imaging/rf02/TypLexMEG/icaanalysis_results/conn_pca/WB_PCA_acrossverts_SL_mippccoh_pc',num2str(which_pc),'.mat'];
% save(out_path,'PCA_SL_ca')
%% over subs
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
%                 [coeff,score,latent,tsquared,explained,mu]=pca(X,'VariableWeights','variance','Centered',true);
%                 exp_mtx_c(vert_no,time_no,band_no,seed_no)=explained(1);
%                 [C_hmuch_c(vert_no,time_no,band_no,seed_no),C_which_c(vert_no,time_no,band_no,seed_no)]=max(corr(X,score(:,1)));
%                 C_c(vert_no,time_no,:,band_no,seed_no)=corr(X,score(:,1));
%                 %abstract
% %                 X=squeeze(Xa(vert_no,time_no,:,:))';
% %                 [coeff,score,latent,tsquared,explained,mu]=pca(X,'VariableWeights','variance','Centered',true);
% %                 exp_mtx_a(vert_no,time_no,band_no,seed_no)=explained(1);
% %                 [C_hmuch_a(vert_no,time_no,band_no,seed_no),C_which_a(vert_no,time_no,band_no,seed_no)]=max(corr(X,score(:,1)));
% %                 C_a(vert_no,time_no,:,band_no,seed_no)=corr(X,score(:,1));
%             end
%         end
%     end
% end
%% over verts- bands
connms={'mi','ppc','coh'};%'mi',
exp_mtx_c=nan(17,2,length(bands),length(seeds));
C_hmuch_c=nan(17,2,length(bands),length(seeds));
C_which_c=nan(17,2,length(bands),length(seeds));
exp_mtx_a=nan(17,2,length(bands),length(seeds));
C_hmuch_a=nan(17,2,length(bands),length(seeds));
C_which_a=nan(17,2,length(bands),length(seeds));
exp_mtx_s=nan(17,2,length(bands),length(seeds));
C_hmuch_s=nan(17,2,length(bands),length(seeds));
C_which_s=nan(17,2,length(bands),length(seeds));
exp_mtx_b=nan(17,2,length(bands),length(seeds));
C_hmuch_b=nan(17,2,length(bands),length(seeds));
C_which_b=nan(17,2,length(bands),length(seeds));
C_c=nan(17,2,length(bands),length(connms),length(seeds));
C_a=nan(17,2,length(bands),length(connms),length(seeds));
C_s=nan(17,2,length(bands),length(connms),length(seeds));
C_b=nan(17,2,length(bands),length(connms),length(seeds));
nbands=length(bands);
nseeds=length(seeds);
ntimes=2;
nsubs=17;
nverts=20484;
PCA_SL_grand=[];
PCA_SL_grand_sb=[];
PCAstd_SL_grand_b=[];
for which_pc=1:3
    which_pc
nconns=3;
for seed_no=1:length(seeds)
    seed=seeds{seed_no};
    seed
    Xc_grand=zeros(nverts,ntimes,nconns,nsubs,nbands);
            Xa_grand=zeros(nverts,ntimes,nconns,nsubs,nbands);
            Xs_grand=zeros(nverts,ntimes,nconns,nsubs,nbands);
            Xb_grand=zeros(2*nverts,ntimes,nconns,nsubs,nbands);
            Xc_grando=zeros(nverts,ntimes,nconns,nsubs,nbands);
            Xa_grando=zeros(nverts,ntimes,nconns,nsubs,nbands);
            Xs_grando=zeros(nverts,ntimes,nconns,nsubs,nbands);
            Xb_grando=zeros(2*nverts,ntimes,nconns,nsubs,nbands);
            scorec_grand=zeros(nverts,ntimes,nconns,nsubs,nbands);
            scorea_grand=zeros(nverts,ntimes,nconns,nsubs,nbands);
            scores_grand=zeros(nverts,ntimes,nconns,nsubs,nbands);
            scoreb_grand=zeros(2*nverts,ntimes,nconns,nsubs,nbands);
    for band_no=1:nbands
        band=bands{band_no};
        band
        load(['/imaging/rf02/TypLexMEG/','3conn_newag_orig_cncrt_',seed,'_',band,'.mat'])
        load(['/imaging/rf02/TypLexMEG/','3conn_newag_orig_abs_',seed,'_',band,'.mat'])
        
        %%%replace outlier in mi
        for xc1=1:size(Xc,2)
                for xc2=1:size(Xc,4)
                    [cc,ccw]=sort(Xc(:,xc1,1,xc2),'descend');
                    [aa,aaw]=sort(Xa(:,xc1,1,xc2),'descend');
                    Xc(ccw(1),xc1,1,xc2)=cc(2);
                    Xa(aaw(1),xc1,1,xc2)=aa(2);
                    
                end
        end
        Xc_grando(:,:,:,:,band_no)=Xc;
        Xa_grando(:,:,:,:,band_no)=Xa;
        Xs_grando(:,:,:,:,band_no)=Xc-Xa;
        Xb_grando(:,:,:,:,band_no)=cat(1,Xc,Xa);
        
        for xcnt1=1:size(Xc,2)
            for xcnt2=1:size(Xc,3)
                for xcnt3=1:size(Xc,4)
                    Xc(:,xcnt1,xcnt2,xcnt3)=boxcox1(Xc(:,xcnt1,xcnt2,xcnt3));
                    Xa(:,xcnt1,xcnt2,xcnt3)=boxcox1(Xa(:,xcnt1,xcnt2,xcnt3));
                end
            end
        end
%         
        Xc_grand(:,:,:,:,band_no)=Xc;
        Xa_grand(:,:,:,:,band_no)=Xa;
        Xs_grand(:,:,:,:,band_no)=Xc-Xa;
        Xb_grand(:,:,:,:,band_no)=cat(1,Xc,Xa);
    end
    for sub_no=1:nsubs
%         sub_no
        for time_no=1:ntimes
            
            for conn_no=1:nbands
                %concrete
                X1=squeeze(Xc_grand(:,time_no,:,sub_no,conn_no));%-Xa(:,time_no,:,sub_no));
                X1o=squeeze(Xc_grando(:,time_no,:,sub_no,conn_no));
                X=X1;%-repmat(mean(X1),size(X1,1),1);
                [coeff,score,latent,tsquared,explained,mu]=pca(X,'VariableWeights','variance','Centered',true);
                exp_mtx_c(sub_no,time_no,conn_no,seed_no)=explained(which_pc);
%                 [C_hmuch_c(sub_no,time_no,conn_no,seed_no),C_which_c(sub_no,time_no,conn_no,seed_no)]=max(corr(X,score(:,which_pc)));
                C_c(sub_no,time_no,conn_no,:,seed_no)=corr(X1,score(:,which_pc));%inv(diag(std(X)))* coeff(:,which_pc);%coeff(:,which_pc);%inv(diag(std(X)))* coeff(:,which_pc);%coeff(:,which_pc);%inv(diag(std(X)))* coeff(:,which_pc);%corr(X,score(:,which_pc));
                scorec_grand(:,time_no,:,sub_no,conn_no)=score;
                %abstract
                X1=squeeze(Xa_grand(:,time_no,:,sub_no,conn_no));
                X1o=squeeze(Xa_grando(:,time_no,:,sub_no,conn_no));
                X=X1;%-repmat(mean(X1),size(X1,1),1);
                [coeff,score,latent,tsquared,explained,mu]=pca(X,'VariableWeights','variance','Centered',true);
                exp_mtx_a(sub_no,time_no,conn_no,seed_no)=explained(which_pc);
%                 [C_hmuch_a(sub_no,time_no,conn_no,seed_no),C_which_a(sub_no,time_no,conn_no,seed_no)]=max(corr(X,score(:,which_pc)));
                C_a(sub_no,time_no,conn_no,:,seed_no)=corr(X1,score(:,which_pc));%inv(diag(std(X)))* coeff(:,which_pc);%coeff(:,which_pc);%inv(diag(std(X)))* coeff(:,which_pc);%coeff(:,which_pc);%inv(diag(std(X)))* coeff(:,which_pc);%corr(X,score(:,which_pc));
                scorea_grand(:,time_no,:,sub_no,conn_no)=score;
                %subtract
                X1=squeeze(Xs_grand(:,time_no,:,sub_no,conn_no));
                X1o=squeeze(Xs_grando(:,time_no,:,sub_no,conn_no));
                X=X1;%-repmat(mean(X1),size(X1,1),1);
                [coeff,score,latent,tsquared,explained,mu]=pca(X,'VariableWeights','variance','Centered',true);
                exp_mtx_s(sub_no,time_no,conn_no,seed_no)=explained(which_pc);
%                 [C_hmuch_s(sub_no,time_no,conn_no,seed_no),C_which_s(sub_no,time_no,conn_no,seed_no)]=max(corr(X,score(:,which_pc)));
                C_s(sub_no,time_no,conn_no,:,seed_no)=corr(X1,score(:,which_pc));%inv(diag(std(X)))* coeff(:,which_pc);%coeff(:,which_pc);%inv(diag(std(X)))* coeff(:,which_pc);%coeff(:,which_pc);%inv(diag(std(X)))* coeff(:,which_pc);%corr(X,score(:,which_pc));
                scores_grand(:,time_no,:,sub_no,conn_no)=score;
                %pooled both
                X1=squeeze(Xb_grand(:,time_no,:,sub_no,conn_no));
                X1o=squeeze(Xb_grando(:,time_no,:,sub_no,conn_no));
                X=X1;%-repmat(mean(X1),size(X1,1),1);
                [coeff,score,latent,tsquared,explained,mu]=pca(X,'VariableWeights','variance','Centered',true);
                exp_mtx_b(sub_no,time_no,conn_no,seed_no)=explained(which_pc);
%                 [C_hmuch_s(sub_no,time_no,conn_no,seed_no),C_which_s(sub_no,time_no,conn_no,seed_no)]=max(corr(X,score(:,which_pc)));
                C_b(sub_no,time_no,conn_no,:,seed_no)=corr(X1,score(:,which_pc));%inv(diag(std(X)))* coeff(:,which_pc);%coeff(:,which_pc);%inv(diag(std(X)))* coeff(:,which_pc);%coeff(:,which_pc);%inv(diag(std(X)))* coeff(:,which_pc);%corr(X,score(:,which_pc));
                scoreb_grand(:,time_no,:,sub_no,conn_no)=score;
            end
        end
    end
end

E=mean(exp_mtx_a,4);
E=squeeze(mean(E,1));
E=squeeze(mean(E,1));

E1=mean(C_a,5);
E1=squeeze(mean(E1,1));
E1=squeeze(mean(E1,1));
PCA_SLa=[E;E1'];

E=mean(exp_mtx_c,4);
E=squeeze(mean(E,1));
E=squeeze(mean(E,1));

E1=mean(C_c,5);
E1=squeeze(mean(E1,1));
E1=squeeze(mean(E1,1));
PCA_SLc=[E;E1'];

E=mean(exp_mtx_s,4);
E=squeeze(mean(E,1));
E=squeeze(mean(E,1));

E1=mean(C_s,5);
E1=squeeze(mean(E1,1));
E1=squeeze(mean(E1,1));
PCA_SLs=[E;E1'];

E=mean(exp_mtx_b,4);
E=squeeze(mean(E,1));
E=squeeze(mean(E,1));
Estd=zeros(size(E));
for estdi=1:length(Estd)
    testd=squeeze(exp_mtx_b(:,:,estdi,:));
    Estd(estdi)=std(testd(:));
end
E1=mean(C_b,5);
E1=squeeze(mean(E1,1));
E1=squeeze(mean(E1,1));
Estd1=zeros(size(E1));
for estdi1=1:size(Estd1,1)
    for estdi2=1:size(Estd1,2)
    testd=squeeze(C_b(:,:,estdi1,estdi2,:));
    Estd1(estdi1,estdi2)=std(testd(:));
    end
end
PCA_SLb=[E;E1'];
PCAstd_SLb=[Estd;Estd1'];

PCA_SL_sb=[PCA_SLb;PCA_SLs];
PCA_SL_ca=[PCA_SLc;PCA_SLa;PCA_SLs];
PCA_SL_grand=[PCA_SL_grand,PCA_SL_ca];
PCA_SL_grand_sb=[PCA_SL_grand_sb,PCA_SL_sb];
PCAstd_SL_grand_b=[PCAstd_SL_grand_b,PCAstd_SLb];

out_path=['/imaging/rf02/TypLexMEG/icaanalysis_results/conn_pca/WB_PCA_acrossverts_SL_allbands_pc',num2str(which_pc),'.mat'];
% save(out_path,'PCA_SL_ca')
end

scoreb_grand1=(scoreb_grand(1:20484,:,:,:,:)+scoreb_grand(20485:end,:,:,:,:))/2;
scoreb_grand1=squeeze(mean(scoreb_grand1,4));
% scoreb_grand1=squeeze(mean(scoreb_grand1,4));
Xb_grand1=(Xb_grand(1:20484,:,:,:,:)+Xb_grand(20485:end,:,:,:,:))/2;
Xb_grand1=squeeze(mean(Xb_grand1,4));
% Xb_grand1=squeeze(mean(Xb_grand1,4));
pc1=squeeze(scoreb_grand1(:,:,1,4));
pc2=squeeze(scoreb_grand1(:,:,2,4));
pc3=squeeze(scoreb_grand1(:,:,3,4));
migb=squeeze(Xb_grand1(:,:,1,4));
cohgb=squeeze(Xb_grand1(:,:,2,4));
ppcgb=squeeze(Xb_grand1(:,:,3,4));
hold on, scatter(pc3(:),migb(:),'y')%figure,
hold on, scatter(pc2(:),migb(:),'r')
hold on, scatter(pc1(:),migb(:),'b')
hold on, scatter(pc3(:),cohgb(:),'y')%figure,
hold on, scatter(pc2(:),cohgb(:),'r')
hold on, scatter(pc1(:),cohgb(:),'b')
hold on, scatter(pc3(:),ppcgb(:),'y')%figure,
hold on, scatter(pc2(:),ppcgb(:),'r')
hold on, scatter(pc1(:),ppcgb(:),'b')

pvmi3=polyfit(pc3(:),migb(:),1);%figure,
pvmi2=polyfit(pc2(:),migb(:),1);%figure,
pvmi1=polyfit(pc1(:),migb(:),1);%figure,

pvcoh3=polyfit(pc3(:),cohgb(:),1);%figure,
pvcoh2=polyfit(pc2(:),cohgb(:),1);%figure,
pvcoh1=polyfit(pc1(:),cohgb(:),1);%figure,

pvppc3=polyfit(pc3(:),ppcgb(:),1);%figure,
pvppc2=polyfit(pc2(:),ppcgb(:),1);%figure,
pvppc1=polyfit(pc1(:),ppcgb(:),1);%figure,

