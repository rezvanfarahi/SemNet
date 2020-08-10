close all
clear
clc
addpath('/home/rf02/rezvan/test1/BCT/2015_01_25 BCT')
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
subjects={'SS1_Meg10_0378',
'SS2_Meg10_0390',
'SS3_Meg11_0026',
'SS4_Meg11_0050',
'SS5_Meg11_0052',
'SS6_Meg11_0069',
'SS9_Meg11_0086',
'SS10_Meg11_0091',
'SS12_Meg11_0096',
'SS13_Meg11_0101',
'SS14_Meg11_0102',
'SS15_Meg11_0112',
'SS16_Meg11_0104',
'SS18_Meg11_0118',
'SS19_Meg11_0131',
'SS20_Meg11_0144',
'SS21_Meg11_0147'
};
% which_parc=1;
% which_conn=2;
outname={'aparc','aparc_mod','aparc2009','aparc2009_mod','RG'};
fname_a={'_ideal_coh.mat','_coh.mat','_ideal_imcoh.mat','_imcoh.mat'};%'aparc_coh_mod.mat';'_coh_mod.mat','_imcoh_mod.mat',
deg_table=struct();
probmat_hubs_all=[];
prob_corr=zeros(5,4);
for which_parc=1:5
    figure,
    for which_conn=1:4
        fname=[label_path,outname{which_parc},fname_a{which_conn}];
        load(fname)
        switch which_conn
            case 2
                wc=1; %con_res_final_ideal;
                
%             case 2
%                 con_res_final=con_mod_final;
            case 1
                con_res_final=con_res_final_ideal;
            case 4
                con_res_final=con_res_final2;
%             case 5
%                 con_res_final=con2_mod_final;
            case 3
                con_res_final=con_res_final_ideal2;
        end
        
        n_subj=size(con_res_final,3);%length(list_all);
        thresh_mat=zeros(size(con_res_final));
        thresh_mat_rand=zeros(size(con_res_final));
        Cof=0.15;
        Deg_a=zeros(n_subj,size(con_res_final,1));
        Deg_rnd=zeros(n_subj,size(con_res_final,1));
        for cnt=1:n_subj
            if sum(which_conn==[1,2])>0
                fname_rnd=[label_path,outname{1},'_',subjects{cnt},'_coh_ideal_max.mat'];
                load(fname_rnd)
                rand_max=coh_ideal_max;
            else
                fname_rnd=[label_path,outname{1},'_',subjects{cnt},'_imcoh_ideal_max.mat'];
                load(fname_rnd)
                rand_max=imcoh_ideal_max;
            end
            sim=sort(rand_max);
            tr_rnd=sim(99);
%             cnt
            con_mat=squeeze(con_res_final(:,:,cnt));
            for cnt1=1:size(con_mat,1)
                for cnt2=1:size(con_mat,2)
                    if cnt2>cnt1
                        con_mat(cnt1,cnt2)=con_mat(cnt2,cnt1);
                    end
                end
            end
            con_mat(abs(con_mat)<tr_rnd)=0;
            conall=con_mat(:);
            [aa,bb]=sort(abs(conall),'descend');
            tr=aa(round(Cof*length(aa)));
            con_mat_rand=-1+2*rand(size(con_mat));
            conallr=con_mat_rand(:);
            [aarnd,bbrnd]=sort(abs(conallr),'descend');
            trrnd=aarnd(round(Cof*length(aarnd)));
            thresh_mat_rand(:,:,cnt)=abs(con_mat_rand)>trrnd;
            thresh_mat_rand(:,:,cnt)=thresh_mat_rand(:,:,cnt)-diag(diag(thresh_mat_rand(:,:,cnt)));
            thresh_mat(:,:,cnt)=abs(con_mat)>tr;
            Deg_a(cnt,:)=sum(thresh_mat(:,:,cnt));
            Deg_rnd(cnt,:)=sum(thresh_mat_rand(:,:,cnt));
        end
        prob_mat=mean(thresh_mat,3);
        prob_mat_rnd=mean(thresh_mat_rand,3);
        prob_tril=tril(prob_mat,-1);
        prob_rnd=sum(prob_tril(:))/((size(prob_tril,1)*(size(prob_tril,1)-1))/2);
        prob_mat(prob_mat<prob_rnd)=0;
%         prob_mat1=double(prob_mat>0);%(3/n_subj)
        deg_vect=mean(Deg_a,1);%sum(prob_mat1);%sum(prob_mat1);
        subplot(2,2,which_conn), hist(deg_vect)
        outmat=[fname(1:end-4),'prob_mat.mat'];
        save(outmat,'prob_mat')
        degvect=[fname(1:end-4),'deg_vect.mat'];
        save(degvect,'deg_vect')
        where_hubs=find(deg_vect>mean(deg_vect+2*std(deg_vect)));
        length(where_hubs)
        where_hubs
        est_deg=deg_vect(where_hubs);
%         est_deg
        est_degr=round(est_deg);
        [pmatsort,pmatsortw]=sort(prob_mat,'descend');
        probmat_hubs=zeros(size(prob_mat));
        hcnt=0;
        for wh=where_hubs
            hcnt=hcnt+1;
            probmat_hubs(pmatsortw(1:est_degr(hcnt),wh),wh)=prob_mat(pmatsortw(1:est_degr(hcnt),wh),wh);
            probmat_hubs(wh,pmatsortw(1:est_degr(hcnt),wh))=prob_mat(pmatsortw(1:est_degr(hcnt),wh),wh);
        end
        outmat=[fname(1:end-4),'prob_mat_hubs.mat'];
        save(outmat,'probmat_hubs')
        if which_conn==1
            ideal_probmat_hubs=probmat_hubs;
            ideal_probmat_hubs(ideal_probmat_hubs<0.5)=0;
        end
        prob_corr(which_parc,which_conn)=corr2(ideal_probmat_hubs,probmat_hubs);
%         probmat_hubs_all(:,:,which_conn)=probmat_hubs;
        
    end
end
close all
%
%     cofi=0;
%     for Cof=node_perc
%         Cof
%         cofi=cofi+1;
%         tr(cnt,cofi)=aa(round(Cof*length(aa)));
%
%         for kk=1:size(plv_mata,3)
%             kk
%             %         pla=squeeze(abs(plv_mata(:,:,kk,:)));
%             %         tra=mean(pla(pla~=0));
%             %
%             %         plc=squeeze(abs(plv_matc(:,:,kk,:)));
%             %         plr=[reshape(plc,1,size(plc,1)*size(plc,2)*size(plc,3)),reshape(pla,1,size(pla,1)*size(pla,2)*size(pla,3))];
%             %         [aa,bb]=sort(plr,'descend');
%             %         tr=aa(round(length(aa)/4));
%             %         tr
%             %         trc=mean(plc(plc~=0));
%             %tr=mean([tra,trc]);
%             for ll=1:size(plv_mata,4)
%
%                 thresh_a(:,:,kk,ll)=threshold_absolute(plv_mata(:,:,kk,ll),tr(cnt,cofi));
%                 deg_a(:,cofi,kk,ll)=degrees_und(thresh_a(:,:,kk,ll));
%
%
%                 %             thresh_s(:,:,kk,ll)=threshold_absolute(plv_mats(:,:,kk,ll),0.5);
%                 %             deg_s(:,kk,ll)=degrees_und(thresh_s(:,:,kk,ll));
%             end
%         end
% % %         deg_cr=zeros(size(deg_c,1),size(deg_c,2),size(deg_c,3),length([6:5:ll-5]));
% % %         deg_ar=zeros(size(deg_a,1),size(deg_a,2),size(deg_a,3),length([6:5:ll-5]));
% % %         cc=0;
% % %         for bc=6:5:ll-5
% % %             cc=cc+1;
% % %             b1=bc-5; b2=b1+10;
% % %             deg_cr(:,:,:,cc)=mean(deg_c(:,:,:,b1:b2),4);
% % %             deg_ar(:,:,:,cc)=mean(deg_a(:,:,:,b1:b2),4);
% % %         end
%         deg_ar=deg_a;
%
%         deg_a1=mean(deg_ar,4);deg_a2=mean(deg_a1,3);
%         Deg_a(cnt,:,cofi)=deg_a2(:,cofi);
%
%     end
%
% %deg_a=squeeze(mean(plv_mata,1));
% %deg_c=squeeze(mean(plv_matc,1));
% %     min(thresh_a(:))
% %     save([in_path,list_all{cnt},'thresholded_25_pca_ppc_abstract.mat'],'thresh_a')
% %     save([in_path,list_all{cnt},'thresholded_25_pca_ppc_concrete.mat'],'thresh_c')
% %     save([in_path,list_all{cnt},'degree_proportional_25_pca_ppc_abstract.mat'],'deg_a')
% %     save([in_path,list_all{cnt},'degree_proportional_25_pca_ppc_concrete.mat'],'deg_c')
% %     save([in_path,list_all{cnt},'degree_mean_pca_ppc_abstract.mat'],'deg_a')
% %     save([in_path,list_all{cnt},'degree_mean_pca_ppc_concrete.mat'],'deg_c')
% %     save([in_path,list_all{cnt},'degree_absthresh_plv_subtract.mat'],'deg_s')
% %
% %
% end
%
% % Degcrat=Deg_c./repmat(mean(Deg_c,2),1,size(Deg_c,2),1);
% % figure,imagesc(node_perc,1:46,squeeze(mean(Degcrat,1))); colorbar
% % degctest=abs(diff(Degcrat,1,3))-0.1;
% % tmap=zeros(size(degctest,2),size(degctest,3));
% % for ii=1:size(tmap,1)
% % for jj=1:size(tmap,2)
% % tmap(ii,jj)=ttest(degctest(:,ii,jj));
% % end
% % end
% % figure, imagesc(node_perc(2:end),1:46,tmap), colormap(gray)
% %
% % Degarat=Deg_a./repmat(mean(Deg_a,2),1,size(Deg_a,2),1);
% % figure,imagesc(node_perc,1:46,squeeze(mean(Degarat,1))); colorbar
% % degatest=abs(diff(Degarat,1,3))-0.1;
% % tmap=zeros(size(degatest,2),size(degatest,3));
% % for ii=1:size(tmap,1)
% % for jj=1:size(tmap,2)
% % tmap(ii,jj)=ttest(degatest(:,ii,jj));
% % end
% % end
% % figure, imagesc(node_perc(2:end),1:46,tmap), colormap(gray)
%
% %Deg_mat=cat(4,Deg_c,Deg_a);Deg_m=mean(Deg_mat,4);
% Deg_m=Deg_a;
% Degarat=Deg_m./repmat(mean(Deg_m,2),1,size(Deg_m,2),1);
% figure,imagesc(node_perc*100,1:size(Degarat,2),squeeze(mean(Degarat,1))); colorbar%,ax1=subplot(1,2,1);colormap(ax1,jet),
% degatest=abs(diff(Degarat,1,3))-0.05;
% tmap=zeros(size(degatest,2),size(degatest,3));
% for ii=1:size(tmap,1)
% for jj=1:size(tmap,2)
% tmap(ii,jj)=ttest(degatest(:,ii,jj));
% end
% end
% figure,imagesc(node_perc(2:end)*100,1:46,tmap),
% colormap(gray)%ax2=subplot(1,2,2);