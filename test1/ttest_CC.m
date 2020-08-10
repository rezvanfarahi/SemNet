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

%fname_a='plv_abstract_cwtmorlet.mat';
%fname_c='plv_concrete_cwtmorlet.mat';
% fname_s='plv_subtract_cwtmorlet.mat';
X_a=zeros(17,64,4,1201);
X_c=zeros(17,64,4,1201);
X_s=zeros(17,64,4,1201);
for cnt=1:length(list_all)
    cnt
    fnamea=[in_path,list_all{cnt},'degree_proportional_25_pca_coh_abstract.mat'];
    load(fnamea)
    
    fnamec=[in_path,list_all{cnt},'degree_proportional_25_pca_coh_concrete.mat'];
    load(fnamec)
    deg_s=deg_c-deg_a;
    X_s(cnt,:,:,:)=deg_s;
    
    deg_a=deg_a-repmat(mean(deg_a(:,:,200:300),3),[1,1,size(deg_a,3)]);
    X_a(cnt,:,:,:)=deg_a;
    deg_c=deg_c-repmat(mean(deg_c(:,:,200:300),3),[1,1,size(deg_c,3)]);
    X_c(cnt,:,:,:)=deg_c;  
    fnamecc=[in_path,list_all{cnt},'degree_contrast_cb_proportional_25_pca_coh.mat'];
    fnameca=[in_path,list_all{cnt},'degree_contrast_ab_proportional_25_pca_coh.mat'];
    fnamecs=[in_path,list_all{cnt},'degree_contrast_ca_proportional_25_pca_coh.mat'];
    save(fnamecc,'deg_c')
    save(fnameca,'deg_a')
    save(fnamecs,'deg_s')

end
deg_a=deg_a(:,:,300:1000);
deg_c=deg_c(:,:,300:1000);
deg_s=deg_s(:,:,300:1000);

% t_a=zeros(size(deg_a));
% t_c=zeros(size(deg_a));
% t_s=zeros(size(deg_a));
% for ii=1:size(X_a,2)
%     ii
%     for jj=1:size(X_a,3)
%         for kk=1:size(X_a,4)
%             [h,p,c,t]=ttest(X_a(:,ii,jj,kk));
%             t_a(ii,jj,kk)=t.tstat;
%             [h,p,c,t]=ttest(X_c(:,ii,jj,kk));
%             t_c(ii,jj,kk)=t.tstat;
%             [h,p,c,t]=ttest(X_s(:,ii,jj,kk));
%             t_s(ii,jj,kk)=t.tstat;
%         end
%     end
% end
t_a_all=size(squeeze(deg_a(1,:,:))); 
t_c_all=size(squeeze(deg_a(1,:,:)));
t_s_all=size(squeeze(deg_a(1,:,:)));
for jj=1:size(X_a,3)
    jj
    for kk=1:size(X_a,4)
        [h,p,c,t]=ttest(mean(X_a(:,:,jj,kk),2));
        t_a_all(jj,kk)=t.tstat;
        [h,p,c,t]=ttest(mean(X_c(:,:,jj,kk),2));
        t_c_all(jj,kk)=t.tstat;
        [h,p,c,t]=ttest(mean(X_s(:,:,jj,kk),2));
        t_s_all(jj,kk)=t.tstat;
    end
end
figure, imagesc(-200:500,1:5,t_s_all), colorbar(), caxis([-max(abs(t_s_all(:))),max(abs(t_s_all(:)))])
xlabel('time (ms)')
ylabel('Gamma                           Beta                               Alpha                           Theta')
% title('t-values: Concrete-Abstract Coherence Average Node Degree')
out_s=['/imaging/rf02/TypLexMEG/icaanalysis_results/jpg/time_resolved/ROI_connectivity/whole_brain_node_degree_coh_concrete_abstract'];
saveas(gcf,[out_s,'.jpg'])
savefig([out_s,'.fig'])

figure, imagesc(-200:500,1:5,t_a_all), colorbar(), caxis([-max(abs(t_a_all(:))),max(abs(t_a_all(:)))])
xlabel('time (ms)')
ylabel('Gamma                           Beta                               Alpha                           Theta')
% title('t-values: Abstract-Baseline Coherence Average Node Degree')
out_s=['/imaging/rf02/TypLexMEG/icaanalysis_results/jpg/time_resolved/ROI_connectivity/whole_brain_node_degree_coh_abstract_baseline'];
saveas(gcf,[out_s,'.jpg'])
savefig([out_s,'.fig'])


figure, imagesc(-200:500,1:5,t_c_all), colorbar(), caxis([-max(abs(t_c_all(:))),max(abs(t_c_all(:)))])
xlabel('time (ms)')
ylabel('Gamma                           Beta                               Alpha                           Theta')
% title('t-values: Concrete-Baseline Coherence Average Node Degree')
out_s=['/imaging/rf02/TypLexMEG/icaanalysis_results/jpg/time_resolved/ROI_connectivity/whole_brain_node_degree_coh_concrete_baseline'];
saveas(gcf,[out_s,'.jpg'])
savefig([out_s,'.fig'])

close all
