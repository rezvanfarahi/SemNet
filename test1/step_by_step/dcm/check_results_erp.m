clear
clc
% close all
addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/')
% addpath('\\cbsu\data\imaging\local\software\spm_cbu_svn\releases\spm12_latest\')

% spm eeg
% spm('defaults','eeg')
cd('/imaging/rf02/TypLexMEG/dcm/latest/')

cwd = '/home/rf02/rezvan/test1/dcm/'

dosubs = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17];
list_all = {'./meg10_0378',
    './meg10_0390',
    './meg11_0026',
    './meg11_0050',
    './meg11_0052',
    './meg11_0069',
    './meg11_0086',
    './meg11_0091',
    './meg11_0096',
    './meg11_0101',
    './meg11_0102',
    './meg11_0112',
    './meg11_0104',
    './meg11_0118',
    './meg11_0131',
    './meg11_0144',
    './meg11_0147',
    };




dcm_path='/imaging/rf02/TypLexMEG/dcm/dipole/';%'/home/rf02/rezvan/test1/step_by_step/dcm/';%'/imaging/rf02/TypLexMEG/dcm/5ROIs_hubvsfull/filtered_250ms_28models_5ROIs/semloc/';%'/imaging/rf02/TypLexMEG/dcm/5ROIs_hubvsfull/250ms_38models_5ROIs/semloc0/';
H_C=[];
H_A=[];
mat_filt=zeros(5,1201,2);
lpfilter = designfilt('lowpassfir', 'FilterOrder',399,'PassbandFrequency', 30, 'StopbandFrequency', 35, 'SampleRate', 1000);
for s = 1:length(dosubs) %parfor
    sess=struct();
    sub    = dosubs(s);
    load(sprintf([dcm_path,list_all{sub}(3:end),'_SemLoc_Evoked_5ROIs_dipole_avg_exttc_50verts_avg.mat']));
%     load(sprintf([dcm_path,list_all{sub}(3:end),'_SemLoc_Evoked_5ROIs_meanCTF_50verts_aifg_forDCM_avg.mat']));
    
    for cc1=1:size(conabs_mat,1)
        for cc2=1:size(conabs_mat,3)
            
        mat_filt(cc1,:,cc2)=(filtfilt(lpfilter,squeeze(conabs_mat(cc1, :, cc2))));
        end
    end
%     mat_filt_ds=zeros(5,240,2);
% jj=0;
% for ii=1:5:1200
% jj=jj+1;
% mat_filt_ds(:,jj,:)=squeeze(mean(mat_filt(:,ii:ii+5,:),2));
% end
    H_C(:,:,s)=squeeze(conabs_mat(:,:,1));
    H_A(:,:,s)=squeeze(conabs_mat(:,:,2));
    
end
asig=mean(H_A,3)';
csig=mean(H_C,3)';
taxis=-500:700;
tmin=400;
tmax=800;
figure,
plot(taxis(tmin:tmax),csig(tmin:tmax,1),'r'), hold on, plot(taxis(tmin:tmax),asig(tmin:tmax,1))
% plot(D(1,500:950,1),'r'), hold on, plot(D(1,500:950,2))
xlabel('time(ms) peri-stimulus')
ylabel('predicted ERP')
title('Anterior Temporal Lobe')

figure,
plot(taxis(tmin:tmax),csig(tmin:tmax,2),'r'), hold on, plot(taxis(tmin:tmax),asig(tmin:tmax,2))
% plot(D(2,500:950,1),'r'), hold on, plot(D(2,500:950,2))
xlabel('time(ms) peri-stimulus')
ylabel('predicted ERP')
title('Inferior Frontal Gyrus')

figure,
plot(taxis(tmin:tmax),csig(tmin:tmax,3),'r'), hold on, plot(taxis(tmin:tmax),asig(tmin:tmax,3))
% plot(D(3,500:950,1),'r'), hold on, plot(D(3,500:950,2))
xlabel('time(ms) peri-stimulus')
ylabel('predicted ERP')
title('Middle Temporal Gyrus')
%
figure,
% plot(D(4,500:750,1),'r'), hold on, plot(D(4,500:750,2))
plot(taxis(tmin:tmax),csig(tmin:tmax,4),'r'), hold on, plot(taxis(tmin:tmax),asig(tmin:tmax,4))
xlabel('time(ms) peri-stimulus')
ylabel('predicted ERP')
title('Angular Gyrus')
%
figure,
plot(taxis(tmin:tmax),csig(tmin:tmax,5),'r'), hold on, plot(taxis(tmin:tmax),asig(tmin:tmax,5))
% plot(D(5,500:950,1),'r'), hold on, plot(D(5,500:950,2))
xlabel('time(ms) peri-stimulus')
ylabel('predicted ERP')
title('Visual Word Form Area')



%        plot(DCM.H{1, 2}(:,end))%H{1,1} concrete, H{1,2} abstract; H{}(:,ii) ii=ROI
%         this_model=struct();
% this_model.fname=sprintf([dcm_path,'DCM_erp_SemLoc_5ROIs_sub%d_mod%d.mat'],sub,n);%DCM.name;