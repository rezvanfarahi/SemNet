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
fname_b='pca_ppc_bothcon_cwtmorlet_resmpled.mat';

fname=[in_path,list_all{1},fname_a];
load(fname)
% fname_s='plv_subtract_cwtmorlet.mat';
n_subj=length(list_all);
node_perc=[0.05:0.05:0.75];

tr=zeros(n_subj,1);
trb=zeros(n_subj,1);

DEG_br=zeros(46,4,2,17);
load('/imaging/rf02/TypLexMEG/icaanalysis_results/jpg/time_resolved/ROI_connectivity/graph/labels_46.mat')
load('/imaging/rf02/TypLexMEG/icaanalysis_results/jpg/time_resolved/ROI_connectivity/graph/node_angles46.mat')
for cnt=1:n_subj
%     cnt
%     fname=[in_path,list_all{cnt},fname_a];
%     load(fname)
%     plv_mata=abs(con_abs(:,:,:,101:end-30));%abs(con_abs(:,:,:,300:end-150));
%     
%     fname=[in_path,list_all{cnt},fname_c];
%     load(fname)
%     plv_matc=abs(con_cncrt(:,:,:,101:end-30));%abs(con_cncrt(:,:,:,300:end-150));
%     
%     
%    
%     conb=cat(5,abs(con_cncrt),abs(con_abs));con_both=mean(conb,5);
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
%     
%     plv_matcr=zeros(size(plv_matc,1),size(plv_matc,2),size(plv_matc,3),2);
%     plv_matar=zeros(size(plv_mata,1),size(plv_mata,2),size(plv_mata,3),2);
%     plv_matbr=zeros(size(plv_matb,1),size(plv_matb,2),size(plv_matb,3),2);
% %     for ii=1:10
% %         plv_matcr(:,:,:,ii)=mean(plv_matc(:,:,:,(ii-1)*10+1:(ii+1)*10),4);
% %         plv_matar(:,:,:,ii)=mean(plv_mata(:,:,:,(ii-1)*10+1:(ii+1)*10),4);
% %         plv_matbr(:,:,:,ii)=mean(plv_matb(:,:,:,(ii-1)*10+1:(ii+1)*10),4);
% %     end
% for ii=1:2
%     plv_matcr(:,:,:,ii)=mean(plv_matc(:,:,:,ii*20+10:(ii+2)*20+10),4);
%     plv_matar(:,:,:,ii)=mean(plv_mata(:,:,:,ii*20+10:(ii+2)*20+10),4);
%     plv_matbr(:,:,:,ii)=mean(plv_matb(:,:,:,ii*20+10:(ii+2)*20+10),4);
% end
% %     for ii=1:10
% %         plv_matcr(:,:,:,ii)=mean(plv_matc(:,:,:,10*ii-5:10*ii+5),4);
% %         plv_matar(:,:,:,ii)=mean(plv_mata(:,:,:,10*ii-5:10*ii+5),4);
% %         plv_matbr(:,:,:,ii)=mean(plv_matb(:,:,:,10*ii-5:10*ii+5),4);
% %     end
% %     for ii=1:5
% %         plv_matcr(:,:,:,ii)=mean(plv_matc(:,:,:,20*ii-10:20*ii+10),4);
% %         plv_matar(:,:,:,ii)=mean(plv_mata(:,:,:,20*ii-10:20*ii+10),4);
% %         plv_matbr(:,:,:,ii)=mean(plv_matb(:,:,:,20*ii-10:20*ii+10),4);
% %     end
%     plv_matb=plv_matbr;
%     plv_matc=plv_matcr;
%     plv_mata=plv_matar;
%     thresh_a=zeros(size(plv_mata));
%     deg_a=zeros(size(plv_mata,2),size(plv_mata,3),size(plv_mata,4));
%     thresh_c=zeros(size(plv_matc));
%     deg_c=zeros(size(plv_matc,2),size(plv_matc,3),size(plv_matc,4));
%     thresh_b=zeros(size(plv_matc));
%     deg_b=zeros(size(plv_matc,2),size(plv_matc,3),size(plv_matc,4));
%     
%     plr=[plv_matc(:);plv_mata(:)];
%     [aa,bb]=sort(plr,'descend');    
%     plrb=plv_matb(:);
%     [aab,bbb]=sort(plrb,'descend'); 
%     
%     Cof=0.15;
%     tr(cnt)=aa(round(Cof*length(aa)));
%     trb(cnt)=aab(round(Cof*length(aab)));
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
%             thresh_b(:,:,kk,ll)=threshold_absolute(plv_matb(:,:,kk,ll),trb(cnt));
%             deg_b(:,kk,ll)=degrees_und(thresh_b(:,:,kk,ll));
%             
%             
%             %             thresh_s(:,:,kk,ll)=threshold_absolute(plv_mats(:,:,kk,ll),0.5);
%             %             deg_s(:,kk,ll)=eigenvector_centrality_und(thresh_s(:,:,kk,ll));
%         end
%     end
%     %deg_mat=cat(4,deg_c,deg_a);deg_b=mean(deg_mat,4);
% %     deg_cr=zeros(size(deg_c,1),size(deg_c,2),length([6:5:size(deg_c,3)-5]));
% %     deg_ar=zeros(size(deg_a,1),size(deg_a,2),length([6:5:size(deg_a,3)-5]));
% %     deg_br=zeros(size(deg_b,1),size(deg_b,2),length([6:5:size(deg_b,3)-5]));
% %     cc=0;
% %     for bc=6:5:size(deg_c,3)-5
% %         cc=cc+1;
% %         b1=bc-5; b2=b1+10;
% %         deg_cr(:,:,cc)=mean(deg_c(:,:,b1:b2),3);
% %         deg_ar(:,:,cc)=mean(deg_a(:,:,b1:b2),3);
% %         deg_br(:,:,cc)=mean(deg_b(:,:,b1:b2),3);
% %     end
%     deg_cr=deg_c;deg_ar=deg_a;deg_br=deg_b;
%     
%     plv_mata1=plv_mata;
%     plv_mata1(thresh_a==0)=0;
%     thresh_ab=thresh_a;
%     thresh_ab(thresh_a>0)=1;
%     strength_a=squeeze(sum(plv_mata1,1));%./squeeze(sum(thresh_ab,1));
%     strength_a(isnan(strength_a))=0;
%     
%     plv_matc1=plv_matc;
%     plv_matc1(thresh_c==0)=0;
%     thresh_cb=thresh_c;
%     thresh_cb(thresh_c>0)=1;
%     strength_c=squeeze(sum(plv_matc1,1));%./squeeze(sum(thresh_cb,1));
%     strength_c(isnan(strength_c))=0;
%     
%     %deg_a=squeeze(mean(plv_mata,1));
%     %deg_c=squeeze(mean(plv_matc,1));
%     min(thresh_a(:))
%     save([in_path,list_all{cnt},fname_b],'con_both')
%     save([in_path,list_all{cnt},'thresholded_15_pca_ppc_abstract_resampled_2win.mat'],'thresh_a')
%     save([in_path,list_all{cnt},'thresholded_15_pca_ppc_concrete_resampled_2win.mat'],'thresh_c')
%     save([in_path,list_all{cnt},'thresholded_15_pca_ppc_both_resampled_2win.mat'],'thresh_b')
%     
%     save([in_path,list_all{cnt},'degrees_proportional_15_pca_ppc_abstract_resampled_2win.mat'],'deg_ar')
%     save([in_path,list_all{cnt},'degrees_proportional_15_pca_ppc_concrete_resampled_2win.mat'],'deg_cr')
%     save([in_path,list_all{cnt},'degrees_proportional_15_pca_ppc_both_resampled_2win.mat'],'deg_br')
%     save([in_path,list_all{cnt},'nodestrength_15_pca_ppc_abstract_resampled_2win.mat'],'strength_a')
%     save([in_path,list_all{cnt},'nodestrength_15_pca_ppc_concrete_resampled_2win.mat'],'strength_c')
%     %         save([in_path,list_all{cnt},'degree_mean_pca_ppc_abstract.mat'],'deg_a')
%     %         save([in_path,list_all{cnt},'degree_mean_pca_ppc_concrete.mat'],'deg_c')
%     %         save([in_path,list_all{cnt},'degree_absthresh_plv_subtract.mat'],'deg_s')
    
    %
    load([in_path,list_all{cnt},'degrees_proportional_15_pca_ppc_concrete_resampled_2win.mat'])
    load([in_path,list_all{cnt},'degrees_proportional_15_pca_ppc_abstract_resampled_2win.mat'])
%     load([in_path,list_all{cnt},'degrees_proportional_15_pca_ppc_both_resampled.mat'])
    degb=cat(4,deg_cr,deg_ar);deg_br=mean(degb,4);
    DEG_br(:,:,:,cnt)=deg_br;
end
deg_brt=DEG_br;%(:,:,40:end,:);
%deg_brr=DEG_br(:,:,1:40,:);

deg_brtz=zeros(size(deg_brt));
%deg_brrz=zeros(size(deg_brr));

% windeg=deg_brt(:,jj,(ii-1)*10+1:(ii+1)*10,kk);
% for ii=1:10
% for jj=1:size(gg,2)
%     for kk=1:size(gg,4)
%         windeg=deg_brt(:,jj,(ii-1)*10+1:(ii+1)*10,kk);
%         deg_brtz(:,jj,ii,kk)=(deg_brt(:,jj,(ii-1)*10+1:(ii+1)*10,kk)-mean(windeg(:)))/std(windeg(:));
%     end
% end
% end
% deg_brtzr=zeros(size(deg_brt,1),size(deg_brt,2),10,size(deg_brt,4));
deg_brtr=deg_brt;%(:,:,2:2:10,:);%zeros(size(deg_brt,1),size(deg_brt,2),10,size(deg_brt,4));
% % % % % % 
% % % % % % deg_brrzp=zeros(size(deg_brtz,1),size(deg_brtz,2),2);
% % % % % % for ii=1:10
% % % % % %     deg_brtr(:,:,ii,:)=mean(deg_brt(:,:,(ii-1)*10+1:(ii+1)*10,:),3);
% % % % % % end
deg_brtzr=zeros(size(deg_brtr));

for jj=1:size(deg_brtr,2) %bands
    for ii=1:size(deg_brtr,3)%time
    jj
    for kk=1:size(deg_brtr,4)%subjs
        windeg=deg_brtr(:,jj,ii,kk);
        deg_brtzr(:,jj,ii,kk)=(deg_brtr(:,jj,ii,kk)-mean(windeg(:)))/std(windeg(:));
%         windeg1=deg_brr(:,jj,:,kk);
%         deg_brrz(:,jj,:,kk)=(deg_brr(:,jj,:,kk)-mean(windeg1(:)))/std(windeg1(:));
    end
    end
end

deg_brtzrm=mean(deg_brtzr,4);
deg_brtzrmz=zeros(size(deg_brtzrm));
for jj=1:size(deg_brtzrm,2)%bands
    for kk=1:size(deg_brtzrm,3)%time
         windeg=deg_brtzrm(:,jj,kk);
         deg_brtzrmz(:,jj,kk)=(deg_brtzrm(:,jj,kk)-mean(windeg(:)))/std(windeg(:));
    end
end
deg_brtzrmzb=deg_brtzrmz>2;
            

deg_brtzrb=deg_brtzr>2;
% deg_brrzb=deg_brrz>1;
deg_brtzrp=mean(deg_brtzrb,4);
% % %  deg_brtzrp=zeros(size(deg_brtzrb,1),size(deg_brtzrb,2),size(deg_brtzrb,3));
% % % % 
% % % % for ii=1:10
% % % %     deg_brtr(:,:,ii,:)=mean(deg_brt(:,:,(ii-1)*10+1:(ii+1)*10,:),3);
% % % %     deg_brtzr(:,:,ii,:)=mean(deg_brtz(:,:,(ii-1)*10+1:(ii+1)*10,:),3);
% % % % end
% % % %deg_brtzr_all=mean(deg_brtzr,4); 
% % % for ii=1:10
% % %     %deg_brtr(:,:,ii,:)=mean(deg_brt(:,:,(ii-1)*10+1:(ii+1)*10,:),3);
% % % for jj=1:size(deg_brtzr,1)%nodes
% % %     for kk=1:size(deg_brtzr,2)%bands
% % %         windeg=deg_brtzrb(jj,kk,ii,:);%mean(deg_brtzb(jj,kk,(ii-1)*10+1:(ii+1)*10,:),3);
% % %         deg_brtzrp(jj,kk,ii)=sum(windeg(:))/length(windeg(:));    
% % %     end
% % % end
% % % end      
%         
% for ii=1:2
% for jj=1:size(deg_brrz,1)%nodes
%     for kk=1:size(deg_brrz,2)%bands
%         windeg=mean(deg_brrzb(jj,kk,(ii-1)*10+1:(ii+1)*10,:),3);
%         deg_brrzp(jj,kk,ii)=sum(windeg(:))/length(windeg(:));    
%     end
% end
% end  
%deg_brtzrp=mean(deg_brtzr,4);
deg_mean=mean(deg_brtr,4);
% deg_final=zeros(size(deg_mean));
% for ii=1:size(deg_final,2)%bands
%     for jj=1:size(deg_final,3)%time
%         
%         deghist=squeeze(deg_brtr(:,ii,2*jj,:));
%         dd=deghist>mean(deghist(:))+3*std(deghist(:));
%         mean(deghist(:))+3*std(deghist(:))
%         dd=squeeze(dd);
%         gg=sum(dd,2)>0;
%         deg_final(gg,ii,jj)=deg_mean(gg,ii,jj);
%     end
% end

deg_final=deg_brtzrp;
deg_p=deg_final;
[aa,bb]=sort(deg_final(:),'descend');
deg_final(deg_final<(mean(aa)+std(aa)))=0;%%(mean(aa)+2*std(aa)))=0;%(mean(aa)+std(aa)))=0;%aa(0.05*length(aa))
%deg_final=deg_brtzrm(:,:,2:2:10);
%deg_final(deg_final<0)=0;
deg_final(deg_final>0)=deg_mean(deg_final>0);
deg_bzmz=deg_brtzrm;
deg_bzmz(deg_final<0)=0;
%deg_final=deg_mean;
%deg_final(deg_p==0)=0;
out_path='/imaging/rf02/TypLexMEG/icaanalysis_results/jpg/time_resolved/ROI_connectivity/graph/ppc_degrees_15_resampled_2win_bothcon2_hubs_prob.mat';
save(out_path,'deg_final')
% out_path='/imaging/rf02/TypLexMEG/icaanalysis_results/jpg/time_resolved/ROI_connectivity/graph/ppc_degrees_zmz_15_resampled_10win_bothcon2_hubs.mat';
% save(out_path,'deg_bzmz')
out_path='/imaging/rf02/TypLexMEG/icaanalysis_results/jpg/time_resolved/ROI_connectivity/graph/ppc_degrees_prob_15_resampled_2win_bothcon2_hubs.mat';
save(out_path,'deg_p')

deg_pca=deg_brtr;
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
