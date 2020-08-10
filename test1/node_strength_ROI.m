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

fname_a='pca_ppc_abstract_cwtmorlet.mat';
fname_c='pca_ppc_concrete_cwtmorlet.mat';
fname=[in_path,list_all{1},fname_a];
load(fname)
% fname_s='plv_subtract_cwtmorlet.mat';
n_subj=length(list_all);
node_perc=[0.05:0.05:0.75];

tr=zeros(n_subj,1);
DEG_br=zeros(46,4,149,17);
load('/imaging/rf02/TypLexMEG/icaanalysis_results/jpg/time_resolved/ROI_connectivity/graph/labels_46.mat')
load('/imaging/rf02/TypLexMEG/icaanalysis_results/jpg/time_resolved/ROI_connectivity/graph/node_angles46.mat')
for cnt=1:n_subj
%     cnt
%     fname=[in_path,list_all{cnt},fname_a];
%     load(fname)
%     plv_mata=abs(con_abs(:,:,:,300:end-150));
%     thresh_a=zeros(size(plv_mata));
%     deg_a=zeros(size(plv_mata,2),size(plv_mata,3),size(plv_mata,4));
%     
%     fname=[in_path,list_all{cnt},fname_c];
%     load(fname)
%     plv_matc=abs(con_cncrt(:,:,:,300:end-150));
%     thresh_c=zeros(size(plv_matc));
%     deg_c=zeros(size(plv_matc,2),size(plv_matc,3),size(plv_matc,4));
%     
%     thresh_b=zeros(size(plv_matc));
%     deg_b=zeros(size(plv_matc,2),size(plv_matc,3),size(plv_matc,4));
%     
%     %     fname=[in_path,list_all{cnt},fname_s];
%     %     load(fname)
%     %     plv_mats=plv_subtract_cwtmorlet;
%     %     thresh_s=zeros(size(plv_mats));
%     %     deg_s=zeros(size(plv_mats,2),size(plv_mats,3),size(plv_mats,4));
%     %
%     for ii=1:size(plv_mata,1)
%         ii
%         for jj=1:size(plv_mata,2)
%             if jj>ii
%                 plv_mata(ii,jj,:,:)=plv_mata(jj,ii,:,:);
%                 plv_matc(ii,jj,:,:)=plv_matc(jj,ii,:,:);
%                 %                 plv_mats(ii,jj,:,:)=plv_mats(jj,ii,:,:);
%             end
%         end
%     end
%     plvb=cat(5,plv_mata,plv_matc);plv_matb=mean(plvb,5);
%     plr=[reshape(plv_matc,1,size(plv_matc,1)*size(plv_matc,2)*size(plv_matc,3)*size(plv_matc,4)),reshape(plv_mata,1,size(plv_mata,1)*size(plv_mata,2)*size(plv_mata,3)*size(plv_mata,4))];
%     [aa,bb]=sort(plr,'descend');
%     Cof=0.15;
%     tr(cnt)=aa(round(Cof*length(aa)));
%     
%     for kk=1:size(plv_mata,3)
%         kk
%         %         pla=squeeze(abs(plv_mata(:,:,kk,:)));
%         %         tra=mean(pla(pla~=0));
%         %
%         %         plc=squeeze(abs(plv_matc(:,:,kk,:)));
%         %         plr=[reshape(plc,1,size(plc,1)*size(plc,2)*size(plc,3)),reshape(pla,1,size(pla,1)*size(pla,2)*size(pla,3))];
%         %         [aa,bb]=sort(plr,'descend');
%         %         tr=aa(round(length(aa)/4));
%         %         tr
%         %         trc=mean(plc(plc~=0));
%         %tr=mean([tra,trc]);
%         for ll=1:size(plv_mata,4)
%             
%             thresh_a(:,:,kk,ll)=threshold_absolute(plv_mata(:,:,kk,ll),tr(cnt));
%             deg_a(:,kk,ll)=degrees_und(thresh_a(:,:,kk,ll));
%             
%             thresh_c(:,:,kk,ll)=threshold_absolute(plv_matc(:,:,kk,ll),tr(cnt));
%             deg_c(:,kk,ll)=degrees_und(thresh_c(:,:,kk,ll));
%             
%             thresh_b(:,:,kk,ll)=threshold_absolute(plv_matb(:,:,kk,ll),tr(cnt));
%             deg_b(:,kk,ll)=degrees_und(thresh_b(:,:,kk,ll));
%             
%             
%             %             thresh_s(:,:,kk,ll)=threshold_absolute(plv_mats(:,:,kk,ll),0.5);
%             %             deg_s(:,kk,ll)=degrees_und(thresh_s(:,:,kk,ll));
%         end
%     end
%     %deg_mat=cat(4,deg_c,deg_a);deg_b=mean(deg_mat,4);
%     deg_cr=zeros(size(deg_c,1),size(deg_c,2),length([6:5:size(deg_c,3)-5]));
%     deg_ar=zeros(size(deg_a,1),size(deg_a,2),length([6:5:size(deg_a,3)-5]));
%     deg_br=zeros(size(deg_b,1),size(deg_b,2),length([6:5:size(deg_b,3)-5]));
%     cc=0;
%     for bc=6:5:size(deg_c,3)-5
%         cc=cc+1;
%         b1=bc-5; b2=b1+10;
%         deg_cr(:,:,cc)=mean(deg_c(:,:,b1:b2),3);
%         deg_ar(:,:,cc)=mean(deg_a(:,:,b1:b2),3);
%         deg_br(:,:,cc)=mean(deg_b(:,:,b1:b2),3);
%     end
%     
%     
%     %deg_a=squeeze(mean(plv_mata,1));
%     %deg_c=squeeze(mean(plv_matc,1));
%     min(thresh_a(:))
%     save([in_path,list_all{cnt},'thresholded_15_pca_ppc_abstract.mat'],'thresh_a')
%     save([in_path,list_all{cnt},'thresholded_15_pca_ppc_concrete.mat'],'thresh_c')
%     save([in_path,list_all{cnt},'thresholded_15_pca_ppc_both.mat'],'thresh_b')
%     
%     save([in_path,list_all{cnt},'degree_proportional_15_pca_ppc_abstract.mat'],'deg_ar')
%     save([in_path,list_all{cnt},'degree_proportional_15_pca_ppc_concrete.mat'],'deg_cr')
%     save([in_path,list_all{cnt},'degree_proportional_15_pca_ppc_both.mat'],'deg_br')
%     %         save([in_path,list_all{cnt},'degree_mean_pca_ppc_abstract.mat'],'deg_a')
%     %         save([in_path,list_all{cnt},'degree_mean_pca_ppc_concrete.mat'],'deg_c')
%     %         save([in_path,list_all{cnt},'degree_absthresh_plv_subtract.mat'],'deg_s')
%     
%     %
    load([in_path,list_all{cnt},'degree_proportional_15_pca_ppc_abstract.mat'])
    load([in_path,list_all{cnt},'degree_proportional_15_pca_ppc_concrete.mat'])
    load([in_path,list_all{cnt},'degree_proportional_15_pca_ppc_both.mat'])
    DEG_br(:,:,:,cnt)=deg_br;
end
deg_brt=DEG_br(:,:,40:end,:);
deg_brr=DEG_br(:,:,1:40,:);

deg_brtz=zeros(size(deg_brt));
deg_brrz=zeros(size(deg_brr));

% windeg=deg_brt(:,jj,(ii-1)*10+1:(ii+1)*10,kk);
% for ii=1:10
% for jj=1:size(gg,2)
%     for kk=1:size(gg,4)
%         windeg=deg_brt(:,jj,(ii-1)*10+1:(ii+1)*10,kk);
%         deg_brtz(:,jj,ii,kk)=(deg_brt(:,jj,(ii-1)*10+1:(ii+1)*10,kk)-mean(windeg(:)))/std(windeg(:));
%     end
% end
% end
deg_brtzr=zeros(size(deg_brt,1),size(deg_brt,2),10,size(deg_brt,4));
deg_brtr=zeros(size(deg_brt,1),size(deg_brt,2),10,size(deg_brt,4));

deg_brrzp=zeros(size(deg_brtz,1),size(deg_brtz,2),2);
for ii=1:10
    deg_brtr(:,:,ii,:)=mean(deg_brt(:,:,(ii-1)*10+1:(ii+1)*10,:),3);
end
for jj=1:size(deg_brtr,2) %bands
    %for ii=1:size(deg_brtr,3)%time
    jj
    for kk=1:size(deg_brtr,4)%subjs
        windeg=deg_brtr(:,jj,:,kk);
        deg_brtzr(:,jj,:,kk)=(deg_brtr(:,jj,:,kk)-mean(windeg(:)))/std(windeg(:));
%         windeg1=deg_brr(:,jj,:,kk);
%         deg_brrz(:,jj,:,kk)=(deg_brr(:,jj,:,kk)-mean(windeg1(:)))/std(windeg1(:));
    %end
    end
end

 deg_brtzrb=deg_brtzr>1;
% deg_brrzb=deg_brrz>1;
 deg_brtzrp=zeros(size(deg_brtzrb,1),size(deg_brtzrb,2),size(deg_brtzrb,3));
% 
% for ii=1:10
%     deg_brtr(:,:,ii,:)=mean(deg_brt(:,:,(ii-1)*10+1:(ii+1)*10,:),3);
%     deg_brtzr(:,:,ii,:)=mean(deg_brtz(:,:,(ii-1)*10+1:(ii+1)*10,:),3);
% end
%deg_brtzr_all=mean(deg_brtzr,4); 
for ii=1:10
    %deg_brtr(:,:,ii,:)=mean(deg_brt(:,:,(ii-1)*10+1:(ii+1)*10,:),3);
for jj=1:size(deg_brtzr,1)%nodes
    for kk=1:size(deg_brtzr,2)%bands
        windeg=deg_brtzrb(jj,kk,ii,:);%mean(deg_brtzb(jj,kk,(ii-1)*10+1:(ii+1)*10,:),3);
        deg_brtzrp(jj,kk,ii)=sum(windeg(:))/length(windeg(:));    
    end
end
end      
%         
% for ii=1:2
% for jj=1:size(deg_brrz,1)%nodes
%     for kk=1:size(deg_brrz,2)%bands
%         windeg=mean(deg_brrzb(jj,kk,(ii-1)*10+1:(ii+1)*10,:),3);
%         deg_brrzp(jj,kk,ii)=sum(windeg(:))/length(windeg(:));    
%     end
% end
% end     
deg_mean=mean(deg_brtr(:,:,2:8,:),4);
deg_final=deg_brtzrp(:,:,2:8);
[aa,bb]=sort(deg_final(:),'descend');
deg_final(deg_final<(mean(aa)+std(aa)))=0;%%(mean(aa)+2*std(aa)))=0;%(mean(aa)+std(aa)))=0;%aa(0.05*length(aa))
deg_final(deg_final>0)=deg_mean(deg_final>0);
out_path='/imaging/rf02/TypLexMEG/icaanalysis_results/jpg/time_resolved/ROI_connectivity/graph/ppc_bothcon_hubs_prob.mat';
save(out_path,'deg_final')
deg_pca=deg_brtr(:,:,2:8,:);
cc=0;
deg_pcaf=[];
for ii=1:size(deg_pca,1)
    for jj=1:size(deg_pca,2)
        for kk=1:size(deg_pca,3)
            if deg_final(ii,jj,kk)>0
                cc=cc+1;
                deg_pcaf(cc,:)=deg_pca(ii,jj,kk,:);
            end
        end
    end
end
[c,s,l]=princomp(deg_pcaf');
cumsum(l)./sum(l);
% out_path='/imaging/rf02/TypLexMEG/icaanalysis_results/jpg/time_resolved/ROI_connectivity/graph/ppc_degforpca_hubs_meanz.mat';
% save(out_path,'deg_pca')
