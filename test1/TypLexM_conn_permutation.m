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

fname_a='pca_coh_abstract_cwtmorlet_resmpled.mat';
fname_c='pca_coh_concrete_cwtmorlet_resmpled.mat';
fname_b='pca_coh_bothcon_cwtmorlet_resmpled.mat';

fname=[in_path,list_all{1},fname_a];
load(fname)
% fname_s='plv_subtract_cwtmorlet.mat';
n_subj=length(list_all);
node_perc=[0.05:0.05:0.75];

tr=zeros(n_subj,1);
trb=zeros(n_subj,1);

CON_sr=zeros(46,46,4,2,17);
load('/imaging/rf02/TypLexMEG/icaanalysis_results/jpg/time_resolved/ROI_connectivity/graph/labels_46.mat')
load('/imaging/rf02/TypLexMEG/icaanalysis_results/jpg/time_resolved/ROI_connectivity/graph/node_angles46.mat')
%%deg_final
load('/imaging/rf02/TypLexMEG/icaanalysis_results/jpg/time_resolved/ROI_connectivity/graph/coh_degrees_15_resampled_2win_bothcon2_hubs_prob.mat')

for cnt=1:n_subj
    cnt
    fname=[in_path,list_all{cnt},fname_a];
    load(fname)
    plv_mata=abs(con_abs(:,:,:,101:end-30));%abs(con_abs(:,:,:,300:end-150));
    
    fname=[in_path,list_all{cnt},fname_c];
    load(fname)
    plv_matc=abs(con_cncrt(:,:,:,101:end-30));%abs(con_cncrt(:,:,:,300:end-150));
    
    
    
    conb=cat(5,abs(con_cncrt),abs(con_abs));con_both=mean(conb,5);
    
    %     fname=[in_path,list_all{cnt},fname_s];
    %     load(fname)
    %     plv_mats=plv_subtract_cwtmorlet;
    %     thresh_s=zeros(size(plv_mats));
    %     deg_s=zeros(size(plv_mats,2),size(plv_mats,3),size(plv_mats,4));
    %
    for ii=1:size(plv_mata,1)
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
    for ii=1:2
        plv_matcr(:,:,:,ii)=mean(plv_matc(:,:,:,ii*20+10:(ii+2)*20+10),4);
        plv_matar(:,:,:,ii)=mean(plv_mata(:,:,:,ii*20+10:(ii+2)*20+10),4);
        plv_matbr(:,:,:,ii)=mean(plv_matb(:,:,:,ii*20+10:(ii+2)*20+10),4);
    end
    CON_sr(:,:,:,:,cnt)=plv_matcr-plv_matar;
end
cnt=0;
dl=deg_final>0;
torig=zeros(size(CON_sr,1),size(CON_sr,2),size(CON_sr,3),size(CON_sr,4));
for ii=1:size(CON_sr,2)%nodes
    ii
    for jj=1:size(CON_sr,3)%bands
        for kk=1:size(CON_sr,4)%time
            if deg_final(ii,jj,kk)>0
                for nn=1:size(CON_sr,2)%nodes
                    cnt=cnt+1;
                    [h1,p1,ci1,s1]=ttest(squeeze(CON_sr(ii,nn,jj,kk,:)));
                    torig(ii,nn,jj,kk)=s1.tstat;
                end
            end
        end
    end
end
torig(isnan(torig))=0;
out_t='/imaging/rf02/TypLexMEG/icaanalysis_results/jpg/time_resolved/ROI_connectivity/graph/uvttest_coh_degrees_15_resampled_2win_subtract.mat';
save(out_t,'torig')


% 
% percheckp=zeros(1,500);
% percheckn=zeros(1,500);
% for percnt=1:500
%     percnt
%     deg_brt=zeros(size(CON_sr));
%     rcoef=sign(randn(1,17));
%     for cntt=1:17
%         if rcoef(cntt)~=0
%             deg_brt(:,:,:,:,cntt)=rcoef(cntt)*CON_sr(:,:,:,:,cntt);%(:,:,40:end,:);
%         else
%             deg_brt(:,:,:,:,cntt)=sign(randn(1,1))*CON_sr(:,:,:,:,cntt);%(:,:,40:end,:);
%         end
%     end
%     
%     tper=zeros(length(sum(dl(:))),1);
%     for ii=1:size(CON_sr,2)%nodes
%         for jj=1:size(CON_sr,3)%bands
%             for kk=1:size(CON_sr,4)%time
%                 if deg_final(ii,jj,kk)>0
%                     for nn=1:size(CON_sr,2)%nodes
%                         cnt=cnt+1;
%                         [h1,p1,ci1,s1]=ttest(squeeze(deg_brt(ii,nn,jj,kk,:)));
%                         tper(cnt)=s1.tstat;
%                     end
%                 end
%             end
%         end
%     end
%     tper(isnan(tper))=0;
%     percheckp(percnt)=max(tper);
%     percheckn(percnt)=min(tper);
%     
% end
% 
% con_sigp=zeros(size(CON_sr,1),size(CON_sr,2),size(CON_sr,3),size(CON_sr,4));
% con_sign=zeros(size(CON_sr,1),size(CON_sr,2),size(CON_sr,3),size(CON_sr,4));
% 
% 
% for ii=1:size(con_sigp,1)
%     ii
%     for jj=1:size(con_sigp,2)
%         for kk=1:size(con_sigp,3)
%             for mm=1:size(con_sigp,4)
%                 con_sigp(ii,jj,kk,mm)=1-sum(torig(ii,jj,kk,mm)>percheckp)/length(torig(ii,jj,kk,mm)>percheckp);
%                 con_sign(ii,jj,kk,mm)=1-sum(torig(ii,jj,kk,mm)<percheckn)/length(torig(ii,jj,kk,mm)<percheckn);
%                 
%             end
%         end
%     end
% end
% posp_log=-10*log(con_sigp);
% negp_log=-10*log(con_sign);