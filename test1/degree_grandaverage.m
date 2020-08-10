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

fname_a='pca_psi_abstract_cwtmorlet.mat';
fname_c='pca_psi_concrete_cwtmorlet.mat';
% fname_s='plv_subtract_cwtmorlet.mat';
A_all=zeros(17,64,64,4,1201);
C_all=zeros(17,64,64,4,1201);

for cnt=1:length(list_all)
    cnt
    fname=[in_path,list_all{cnt},fname_a];
    load(fname)
    plv_mata=con_abs;
    thresh_a=zeros(size(plv_mata));
    deg_a=zeros(size(plv_mata,2),size(plv_mata,3),size(plv_mata,4));
    
    fname=[in_path,list_all{cnt},fname_c];
    load(fname)
    plv_matc=con_cncrt;
    thresh_c=zeros(size(plv_matc));
    deg_c=zeros(size(plv_matc,2),size(plv_matc,3),size(plv_matc,4));
    
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
    A_all(cnt,:,:,:,:)=plv_mata;
    C_all(cnt,:,:,:,:)=plv_matc;
    
%     for kk=1:size(plv_mata,3)
%         kk
%         pla=squeeze(plv_mata(:,:,kk,:));
%         tra=mean(pla(pla~=0));
%         
%         plc=squeeze(plv_matc(:,:,kk,:));
%         trc=mean(plc(plc~=0));
%         tr=mean([tra,trc]);
%         for ll=1:size(plv_mata,4)
%             
%             thresh_a(:,:,kk,ll)=threshold_absolute(plv_mata(:,:,kk,ll),tra);
%             deg_a(:,kk,ll)=degrees_und(thresh_a(:,:,kk,ll));
%            
%             thresh_c(:,:,kk,ll)=threshold_absolute(plv_matc(:,:,kk,ll),trc);
%             deg_c(:,kk,ll)=degrees_und(thresh_c(:,:,kk,ll));
%             
% %             thresh_s(:,:,kk,ll)=threshold_absolute(plv_mats(:,:,kk,ll),0.5);
% %             deg_s(:,kk,ll)=degrees_und(thresh_s(:,:,kk,ll));
%         end
%     end
%     save([in_path,list_all{cnt},'degree_absolute_pca_ppc_abstract.mat'],'deg_a')
%     save([in_path,list_all{cnt},'degree_absolute_pca_ppc_concrete.mat'],'deg_c')
% %     save([in_path,list_all{cnt},'degree_absthresh_plv_subtract.mat'],'deg_s')
    
    
end
