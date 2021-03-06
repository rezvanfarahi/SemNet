clear
clc
close all
addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/')
% addpath('\\cbsu\data\imaging\local\software\spm_cbu_svn\releases\spm12_latest\')

% spm eeg
% spm('defaults','eeg')
cd('/imaging/rf02/TypLexMEG/dcm/5ROIs_hubvsfull/')

cwd = '/home/rf02/rezvan/test1/dcm'

dosubs = [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19];
% list_all = {'./meg10_0378',
%     './meg10_0390',
%     './meg11_0026',
%     './meg11_0050',
%     './meg11_0052',
%     './meg11_0069',
%     './meg11_0086',
%     './meg11_0091',
%     './meg11_0096',
%     './meg11_0101',
%     './meg11_0102',
%     './meg11_0112',
%     './meg11_0104',
%     './meg11_0118',
%     './meg11_0131',
%     './meg11_0144',
%     './meg11_0147',
%     };




dcm_path='/imaging/rf02/Semnet/semnet4semloc/dcm/maxCTF/LD/oldreg_filtered_450ms_28models_5ROIs/maxCTF_ERP_dtr0/';%'/imaging/rf02/TypLexMEG/dcm/latest/150ms_57models_5ROIs/semloc/';%semloc/avg_allverts_dtr0_dip/ctf/';%'/imaging/rf02/TypLexMEG/dcm/5ROIs_hubvsfull/250ms_38models_5ROIs/semloc0/';
outpath_anova='/imaging/rf02/Semnet/semnet4semloc/dcm/SemNet_LD_predicted_450ms_AG.mat';
%'/imaging/rf02/TypLexMEG/dcm/5ROIs_hubvsfull/filtered_250ms_28models_5ROIs/semloc/maxCTF_ERP_dtr0/simulated/inverted/model3/';%
H_C=[];
H_A=[];
tmax=451;

for s = 1:length(dosubs) %parfor
    sess=struct();
    sub    = dosubs(s);
    
   ncnt=0;
    for n=[14,16,45:46]%[2,4,39:40]%[30,32]%[2,4,39:40,6,8,41:42,10,12,43:44,14,16,45:46]%,18:2:36,37:38]%[2:2:36,37:38]%1:38%numel(Model)
        ncnt=ncnt+1;
        sub
        n

        load(sprintf([dcm_path,'DCM_erpf_ConEmot_LD_oldreg_5ROIs_avg_maxCTF_sub%d_mod%d.mat'],sub,n));%DCM_erpf_SemLoc_5ROIs_mpowavg_sub%d_mod%d.mat
%         DCM_erpf_SemLoc_5ROIs_avg_allverts_pifg_sub11_mod8
        H_C(:,:,ncnt,s)=DCM.H{1,1};%DCM.xY.y{1,1};%
        H_A(:,:,ncnt,s)=DCM.H{1,2};%DCM.xY.y{1,2};%
    end
end
H_CA=zeros(size(H_C,1),size(H_C,2),size(H_C,4),2);
H_CA(:,:,:,1)=mean(H_C,3);
H_CA(:,:,:,2)=mean(H_A,3);
save(outpath_anova,'H_CA')

asig=mean(H_A(1:tmax,:,:,:),4);%H_A(:,:,4,7);%
asig=mean(asig,3);%squeeze(asig(:,:,4));%
csig=mean(H_C(1:tmax,:,:,:),4);%H_C(:,:,4,7);%
csig=mean(csig,3);%squeeze(csig(:,:,4));%

figure,
plot(csig(:,1),'r','LineWidth',3), hold on, plot(asig(:,1),'LineWidth',3)
% plot(D(1,500:950,1),'r'), hold on, plot(D(1,500:950,2))
xlabel('time(ms) peri-stimulus')
ylabel('predicted ERP')
title('Anterior Temporal Lobe')
% legend('Concrete','Abstract')

figure,
plot(csig(:,2),'r','LineWidth',3), hold on, plot(asig(:,2),'LineWidth',3)
% plot(D(2,500:950,1),'r'), hold on, plot(D(2,500:950,2))
xlabel('time(ms) peri-stimulus')
ylabel('predicted ERP')
title('Inferior Frontal Gyrus')
% legend('Concrete','Abstract')


figure,
plot(csig(:,3),'r','LineWidth',3), hold on, plot(asig(:,3),'LineWidth',3)
% plot(D(3,500:950,1),'r'), hold on, plot(D(3,500:950,2))
xlabel('time(ms) peri-stimulus')
ylabel('predicted ERP')
title('Middle Temporal Gyrus')
% legend('Concrete','Abstract')

% 
figure,
% plot(D(4,500:750,1),'r'), hold on, plot(D(4,500:750,2))
plot(csig(:,4),'r','LineWidth',3), hold on, plot(asig(:,4),'LineWidth',3)
xlabel('time(ms) peri-stimulus')
ylabel('predicted ERP')
title('Angular Gyrus')
% legend('Concrete','Abstract')

% 
figure,
plot(csig(:,5),'r','LineWidth',3), hold on, plot(asig(:,5),'LineWidth',3)
% plot(D(5,500:950,1),'r'), hold on, plot(D(5,500:950,2))
xlabel('time(ms) peri-stimulus')
ylabel('predicted ERP')
title('Visual Word Form Area')
% legend('Concrete','Abstract')


cd('/imaging/rf02/Semnet/semnet4semloc/dcm/')

%        plot(DCM.H{1, 2}(:,end))%H{1,1} concrete, H{1,2} abstract; H{}(:,ii) ii=ROI
%         this_model=struct();
% this_model.fname=sprintf([dcm_path,'DCM_erp_SemLoc_5ROIs_sub%d_mod%d.mat'],sub,n);%DCM.name;