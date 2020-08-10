close all
clear
clc
addpath('/home/rf02/rezvan/test1/BCT/2015_01_25 BCT')
in_path='/imaging/rf02/TypLexMEG/';

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

fname_a='pca_ppc_abstract_cwtmorlet_resmpled.mat';
fname_c='pca_ppc_concrete_cwtmorlet_resmpled.mat';
fname=[in_path,list_all{1},fname_a];
load(fname)
% fname_s='plv_subtract_cwtmorlet.mat';
n_subj=length(list_all);
node_perc=[0.05:0.05:0.75];
Deg_c=zeros(n_subj,size(con_abs,1),length(node_perc));
Deg_a=zeros(n_subj,size(con_abs,1),length(node_perc));
Deg_b=zeros(n_subj,size(con_abs,1),length(node_perc));

tr=zeros(n_subj,length(node_perc));
trb=zeros(n_subj,length(node_perc));


for cnt=1:n_subj
   cnt
    fname=[in_path,list_all{cnt},fname_a];
    load(fname)
    plv_mata=abs(con_abs(:,:,:,100:end));%abs(con_abs(:,:,:,300:end-150));
    
    fname=[in_path,list_all{cnt},fname_c];
    load(fname)
    plv_matc=abs(con_cncrt(:,:,:,100:end));%abs(con_cncrt(:,:,:,300:end-150));
    
   
    conb=cat(5,abs(con_cncrt),abs(con_abs));con_both=mean(conb,5);
    
    %     fname=[in_path,list_all{cnt},fname_s];
    %     load(fname)
    %     plv_mats=plv_subtract_cwtmorlet;
    %     thresh_s=zeros(size(plv_mats));
    %     deg_s=zeros(size(plv_mats,2),size(plv_mats,3),size(plv_mats,4));
    %
    for ii=1:size(plv_mata,1)
        ii
        for jj=1:size(plv_mata,2)
            if jj>ii
                plv_mata(ii,jj,:,:)=plv_mata(jj,ii,:,:);
                plv_matc(ii,jj,:,:)=plv_matc(jj,ii,:,:);
                %                 plv_mats(ii,jj,:,:)=plv_mats(jj,ii,:,:);
            end
        end
    end
    plvb=cat(5,plv_mata,plv_matc);plv_matb=mean(plvb,5);
    
    plv_matcr=zeros(size(plv_matc,1),size(plv_matc,2),size(plv_matc,3),2);
    plv_matar=zeros(size(plv_mata,1),size(plv_mata,2),size(plv_mata,3),2);
    plv_matbr=zeros(size(plv_matb,1),size(plv_matb,2),size(plv_matb,3),2);
%     for ii=1:10
%         plv_matcr(:,:,:,ii)=mean(plv_matc(:,:,:,(ii-1)*10+1:(ii+1)*10),4);
%         plv_matar(:,:,:,ii)=mean(plv_mata(:,:,:,(ii-1)*10+1:(ii+1)*10),4);
%         plv_matbr(:,:,:,ii)=mean(plv_matb(:,:,:,(ii-1)*10+1:(ii+1)*10),4);
%     end
% for ii=1:10
%         plv_matcr(:,:,:,ii)=mean(plv_matc(:,:,:,10*ii-5:10*ii+5),4);
%         plv_matar(:,:,:,ii)=mean(plv_mata(:,:,:,10*ii-5:10*ii+5),4);
%         plv_matbr(:,:,:,ii)=mean(plv_matb(:,:,:,10*ii-5:10*ii+5),4);
%     end
for ii=1:2
    plv_matcr(:,:,:,ii)=mean(plv_matc(:,:,:,ii*20+10:(ii+2)*20+10),4);
    plv_matar(:,:,:,ii)=mean(plv_mata(:,:,:,ii*20+10:(ii+2)*20+10),4);
    plv_matbr(:,:,:,ii)=mean(plv_matb(:,:,:,ii*20+10:(ii+2)*20+10),4);
end
    plv_matb=plv_matbr;
    plv_matc=plv_matcr;
    plv_mata=plv_matar;
    thresh_a=zeros(size(plv_mata));
    deg_a=zeros(size(plv_mata,2),size(plv_mata,3),size(plv_mata,4));
    thresh_c=zeros(size(plv_matc));
    deg_c=zeros(size(plv_matc,2),size(plv_matc,3),size(plv_matc,4));
    thresh_b=zeros(size(plv_matb));
    deg_b=zeros(size(plv_matb,2),size(plv_matb,3),size(plv_matb,4));
    
    plr=[plv_matc(:);plv_mata(:)];
    [aa,bb]=sort(plr,'descend');    
    plrb=plv_matb(:);
    [aab,bbb]=sort(plrb,'descend'); 
    
    cofi=0;
    for Cof=node_perc
        Cof
        cofi=cofi+1;
        tr(cnt,cofi)=aa(round(Cof*length(aa)));
        trb(cnt,cofi)=aab(round(Cof*length(aab)));
    
        for kk=1:size(plv_mata,3)
            kk
            %         pla=squeeze(abs(plv_mata(:,:,kk,:)));
            %         tra=mean(pla(pla~=0));
            %
            %         plc=squeeze(abs(plv_matc(:,:,kk,:)));
            %         plr=[reshape(plc,1,size(plc,1)*size(plc,2)*size(plc,3)),reshape(pla,1,size(pla,1)*size(pla,2)*size(pla,3))];
            %         [aa,bb]=sort(plr,'descend');
            %         tr=aa(round(length(aa)/4));
            %         tr
            %         trc=mean(plc(plc~=0));
            %tr=mean([tra,trc]);
            for ll=1:size(plv_mata,4)
                
                thresh_a(:,:,kk,ll)=threshold_absolute(plv_mata(:,:,kk,ll),tr(cnt,cofi));
                deg_a(:,cofi,kk,ll)=degrees_und(thresh_a(:,:,kk,ll));
                
                thresh_c(:,:,kk,ll)=threshold_absolute(plv_matc(:,:,kk,ll),tr(cnt,cofi));
                deg_c(:,cofi,kk,ll)=degrees_und(thresh_c(:,:,kk,ll));
                
                thresh_b(:,:,kk,ll)=threshold_absolute(plv_matb(:,:,kk,ll),trb(cnt,cofi));
                deg_b(:,cofi,kk,ll)=degrees_und(thresh_b(:,:,kk,ll));
                
                %             thresh_s(:,:,kk,ll)=threshold_absolute(plv_mats(:,:,kk,ll),0.5);
                %             deg_s(:,kk,ll)=degrees_und(thresh_s(:,:,kk,ll));
            end
        end
% %         deg_cr=zeros(size(deg_c,1),size(deg_c,2),size(deg_c,3),length([6:5:ll-5]));
% %         deg_ar=zeros(size(deg_a,1),size(deg_a,2),size(deg_a,3),length([6:5:ll-5]));
% %         cc=0;
% %         for bc=6:5:ll-5
% %             cc=cc+1;
% %             b1=bc-5; b2=b1+10;
% %             deg_cr(:,:,:,cc)=mean(deg_c(:,:,:,b1:b2),4);
% %             deg_ar(:,:,:,cc)=mean(deg_a(:,:,:,b1:b2),4);
% %         end
        deg_cr=deg_c;
        deg_ar=deg_a;
        deg_br=deg_b;
        deg_c1=mean(deg_cr,4);deg_c2=mean(deg_c1,3);
        Deg_c(cnt,:,cofi)=deg_c2(:,cofi);
        deg_a1=mean(deg_ar,4);deg_a2=mean(deg_a1,3);
        Deg_a(cnt,:,cofi)=deg_a2(:,cofi);
        deg_b1=mean(deg_br,4);deg_b2=mean(deg_b1,3);
        Deg_b(cnt,:,cofi)=deg_b2(:,cofi);
    end
    
%deg_a=squeeze(mean(plv_mata,1));
%deg_c=squeeze(mean(plv_matc,1));
%     min(thresh_a(:))
%     save([in_path,list_all{cnt},'thresholded_25_pca_ppc_abstract.mat'],'thresh_a')
%     save([in_path,list_all{cnt},'thresholded_25_pca_ppc_concrete.mat'],'thresh_c')
%     save([in_path,list_all{cnt},'degree_proportional_25_pca_ppc_abstract.mat'],'deg_a')
%     save([in_path,list_all{cnt},'degree_proportional_25_pca_ppc_concrete.mat'],'deg_c')
%     save([in_path,list_all{cnt},'degree_mean_pca_ppc_abstract.mat'],'deg_a')
%     save([in_path,list_all{cnt},'degree_mean_pca_ppc_concrete.mat'],'deg_c')
%     save([in_path,list_all{cnt},'degree_absthresh_plv_subtract.mat'],'deg_s')
%
%
end

% Degcrat=Deg_c./repmat(mean(Deg_c,2),1,size(Deg_c,2),1);
% figure,imagesc(node_perc,1:46,squeeze(mean(Degcrat,1))); colorbar
% degctest=abs(diff(Degcrat,1,3))-0.1;
% tmap=zeros(size(degctest,2),size(degctest,3));
% for ii=1:size(tmap,1)
% for jj=1:size(tmap,2)
% tmap(ii,jj)=ttest(degctest(:,ii,jj));
% end
% end
% figure, imagesc(node_perc(2:end),1:46,tmap), colormap(gray)
% 
% Degarat=Deg_a./repmat(mean(Deg_a,2),1,size(Deg_a,2),1);
% figure,imagesc(node_perc,1:46,squeeze(mean(Degarat,1))); colorbar
% degatest=abs(diff(Degarat,1,3))-0.1;
% tmap=zeros(size(degatest,2),size(degatest,3));
% for ii=1:size(tmap,1)
% for jj=1:size(tmap,2)
% tmap(ii,jj)=ttest(degatest(:,ii,jj));
% end
% end
% figure, imagesc(node_perc(2:end),1:46,tmap), colormap(gray)

%Deg_mat=cat(4,Deg_c,Deg_a);Deg_m=mean(Deg_mat,4);
Deg_m=Deg_b;
Degarat=Deg_m./repmat(mean(Deg_m,2),1,size(Deg_m,2),1);
figure,imagesc(node_perc*100,1:46,squeeze(mean(Degarat,1))); colorbar%,ax1=subplot(1,2,1);colormap(ax1,jet),
degatest=abs(diff(Degarat,1,3))-0.05;
tmap=zeros(size(degatest,2),size(degatest,3));
for ii=1:size(tmap,1)
for jj=1:size(tmap,2)
tmap(ii,jj)=ttest(degatest(:,ii,jj));
end
end
figure,imagesc(node_perc(2:end)*100,1:46,tmap), colormap(gray)%ax2=subplot(1,2,2);