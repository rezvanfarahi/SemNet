% demo for creating an SPM M/EEG dataset from arbitrary data using
% conversion of simple Fieldtrip raw data struct.
% SPM8 internal format is quite complex but is transparent to the user
% as meeg object's methods take care of maintaining its consistency.
% The most straightforward way to convert arbitrary data that is available
% as an ASCII file or *.mat file with some variables to SPM8 is to create
% a quite simple Fieldtrip raw data struct and then use SPM's
% spm_eeg_ft2spm to convert this struct to SPM8 file. Missing information
% can then be supplemented using meeg methods and SPM functions.
% Fieldtrip raw struct must contain the following fields:
% .fsample - sampling rate (Hz)
% .trial - cell array of matrices with identical dimensions channels x time
% .time - cell array of time vectors (in sec), the same length as the
%         second dimension of the data. For SPM8 they must be identical.
% .label - cell array of strings, list of channel labels. Same length as
%         the first dimension of the data.
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_convert_arbitrary_data.m 5404 2013-04-12 15:08:57Z vladimir $
clear
clc
close all
% Initialize SPM
%--------------------------------------------------------------------------
addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r6906/')
% spm eeg
% spm('defaults','EEG');


% Some details about the data
%--------------------------------------------------------------------------
q=1;
Nchannels = 5;
Nsamples  = round(1000/q)+1;
Ntrials   = 2;
TimeOnset = -0.4; % in sec
Fsample = 1000/q;

% chlabels = {
%     
% 'lSMG'
% 'lIFG'
% 'lMTG'
% 'lATL'
% 'lWFA'
% 
% };
chlabels = {
    
'lATL'
'lIFG'
'lMTG'
'lAG'
'lWFA'

};
input_path='/imaging/rf02/TypLexMEG/dcm/maxCTF/';%'/imaging/rf02/TypLexMEG/dcm/centre_of_mass/';%'/imaging/rf02/TypLexMEG/dcm/centre_of_mass/';%'/imaging/rf02/TypLexMEG/dcm/maxCTF/';%'/imaging/rf02/TypLexMEG/dcm/dipole/';
outpath_anova='/imaging/rf02/Semnet/semnet4semloc/dcm/SemLoc_SD_500wins_unsigned.mat';

list_all = {'/meg10_0378',
    '/meg10_0390',
    '/meg11_0026',
    '/meg11_0050',
    '/meg11_0052',
    '/meg11_0069',
    '/meg11_0086',
    '/meg11_0091',
    '/meg11_0096',
    '/meg11_0101',
    '/meg11_0102',
    '/meg11_0112',
    '/meg11_0104',
    '/meg11_0118',
    '/meg11_0131',
    '/meg11_0144',
    '/meg11_0147',
    };
% define the output file name
Data=[];
Dataf=[];
lpfilter = designfilt('lowpassfir', 'FilterOrder',399,'PassbandFrequency', 30, 'StopbandFrequency', 35, 'SampleRate', 1000);

for cnt=1:length(list_all)
    %--------------------------------------------------------------------------
    cnt
    
    % create data array
    %--------------------------------------------------------------------------
    fname_a='_SemLoc_Evoked_unsigned_5ROIs_mean_maxCTF_forDCM_avg';%'_TypLex_SL_Evoked_5ROIs_meanflip_maxCTF_forDCM_avg';%'_SemLoc_Evoked_5ROIs_maxCTF_mean_forDCM_avg';%'_SemLoc_Evoked_5ROIs_meanflip_maxCTF_forDCM_avg_aifg';%'_SemLoc_Evoked_5ROIs_cmass_forDCM';%'_SemLoc_Evoked_5ROIs_cmass_forDCM';%'_SemLoc_Evoked_5ROIs_maxCTF_forDCM_avg';%'_SemLoc_Evoked_5ROIs_relctfdip_avgnoflip_exttc_allverts_avg';%'_SemLoc_Evoked_5ROIs_dipole_pifg_avg_snr3_exttc_allverts_avg';%'_SemLoc_Evoked_5ROIs_meanCTF_50verts_aifg_forDCM_avg';%_SemLoc_Evoked_5ROIs_meanCTF_forDCM_avg.mat
                                                                   
    fname_in=[input_path,list_all{cnt},fname_a,'.mat'];
    this_dir=input_path;%'/imaging/rf02/TypLexMEG/dcm/5ROIs_hubvsfull/avg_50verts_dipole/';
    if ~exist(this_dir,'dir')
        mkdir(this_dir)
    end
%     fname_out = [this_dir,list_all{cnt},fname_a,'_filtered_cnvrt.mat'];
    load(fname_in)
    data = conabs_mat;%-repmat((mean(conabs_mat(:,:,:),2)),1,size(conabs_mat,2),1);%randn([NROIs, Ntimes, Nconditions]);
    Data(:,:,:,cnt)=data;
    % create the time axis (should be the same for all trials)
    %--------------------------------------------------------------------------
    timeaxis = [0:(Nsamples-1)]./Fsample + TimeOnset;
    
    % Create the Fieldtrip raw struct
    %--------------------------------------------------------------------------
    
    ftdata = [];
    
    for i = 1:Ntrials
        %     if data(1,650,i)<0
        %         data(1,:,i)=-data(1,:,i);
        %     end
        %     if data(2,650,i)<0
        %         data(2,:,i)=-data(2,:,i);
        %     end
        thisdata_filtered=[];
        for roicnt=1:size(data,1)
%             q=5;
        thisdata_filtered(roicnt,:)=resample(filtfilt(lpfilter,squeeze(data(roicnt, :, i))),1,q);%squeeze(data(roicnt, :, i));%filtfilt(lpfilter,squeeze(data(roicnt, :, i)));
        end
        Dataf(:,:,i,cnt)=thisdata_filtered(:,100/q:1100/q);
        ftdata.trial{i} = thisdata_filtered(:,100/q:1100/q);%squeeze(data(:, :, i));%thisdata_filtered;%squeeze(data(:, :, i));%*1e10
        ftdata.time{i} = timeaxis;
    end
    
    
    ftdata.fsample = Fsample;
    ftdata.label = chlabels;
    ftdata.label = ftdata.label(:);
    
    % Convert the ftdata struct to SPM M\EEG dataset
    %--------------------------------------------------------------------------
%     D = spm_eeg_ft2spm(ftdata, fname_out);
%     
%     % Examples of providing additional information in a script
%     % [] comes instead of an index vector and means that the command
%     % applies to all channels/all trials.
%     %--------------------------------------------------------------------------
%     D = type(D, 'single');                        % Sets the dataset type
%     D = chantype(D, ':', 'LFP');                   % Sets the channel type
%     D = conditions(D, 1, 'Conctere');  % Sets the condition label
%     D = conditions(D, 2, 'Abstract');
%     
%     % save
%     %--------------------------------------------------------------------------
%     save(D);
end

% Datafiltered
% lpfilter = designfilt('lowpassfir', 'FilterOrder',399,'PassbandFrequency', 30, 'StopbandFrequency', 35, 'SampleRate', 1000);
% ddf=filtfilt(lpfilter,dd);
D1=Dataf(:,:,:,7);%mean(Dataf,4);%abs(squeeze(Dataf(:,:,:,8)));%mean(Dataf,4);
Dstd=zeros(size(Dataf,1),size(Dataf,2),size(Dataf,3));
for dcnt1=1:size(Dstd,1)
    for dcnt2=1:size(Dstd,2)
        for dcnt3=1:size(Dstd,3)
            Dstd(dcnt1,dcnt2,dcnt3)=std(squeeze(Data(dcnt1,dcnt2,dcnt3,:)));
        end
    end
end
Dp=zeros(size(Dataf,1),size(Dataf,2));
Dt=zeros(size(Dataf,1),size(Dataf,2));

for dcnt1=1:size(Dstd,1)
    for dcnt2=1:size(Dstd,2)
            [h,Dp(dcnt1,dcnt2),ci,tstat1]=ttest(squeeze(Dataf(dcnt1,dcnt2,2,:)));
            Dt(dcnt1,dcnt2)=tstat1.tstat;
    end
end
% ntimewins=4;
% Dpm=zeros(size(Dataf,1),ntimewins);
% for ii=1:ntimewins
%     Datafm=squeeze(mean((Dataf(:,(300+50*ii)/q:(350+50*ii)/q,:,:)),2));%5x2x17
%     for dcnt1=1:size(Dpm,1)
%             [h,Dpm(dcnt1,ii)]=ttest(squeeze(Datafm(dcnt1,1,:))-squeeze(Datafm(dcnt1,2,:)));
%     end
% end

ntimewins=500;
Dpm=zeros(size(Dataf,1),ntimewins);
Datafm_wins=zeros(size(Dataf,1),ntimewins,size(Dataf,3),size(Dataf,4));
for ii=1:ntimewins
    Datafm=squeeze(mean((Dataf(:,(350+1*(ii-1))/q:(350+1*ii)/q,:,:)),2));%5x2x17
    Datafm_wins(:,ii,:,:)=Datafm;
    for dcnt1=1:size(Dpm,1)
            [h,Dpm(dcnt1,ii)]=ttest(squeeze(Datafm(dcnt1,1,:))-squeeze(Datafm(dcnt1,2,:)));
    end
end
save(outpath_anova,'Datafm_wins')


Dse=Dstd/sqrt(size(Dataf,4));
ts=tinv([0.025,0.975],size(Dataf,4)-1);
Dci=Dse*ts(2);
taxis=-400:q:600;
tmin=350/q;
tmax=850/q;
% Dstd=std(Dataf,4);
figure,   
plot(taxis(tmin:tmax),D1(1,tmin:tmax,1),'r','LineWidth',3), hold on, plot(taxis(tmin:tmax),D1(1,tmin:tmax,2),'LineWidth',3)
% plot(Dci(1,500:750,1)+D(1,500:750,1),'--r'),plot(Dci(1,500:750,2)+D(1,500:750,2),'--')
% plot(-Dci(1,500:750,1)+D(1,500:750,1),'--r'),plot(-Dci(1,500:750,2)+D(1,500:750,2),'--')
% h=fill([1:251,1:251],[D(1,500:750,1),Dci(1,500:750,1)+D(1,500:750,1)],'r');
% h=area([D(1,500:750,1);Dci(1,500:750,1)+D(1,500:750,1)]);
% set(h,'facealpha',0.5)
% plot(D(1,500:950,1),'r'), hold on, plot(D(1,500:950,2))
xlabel('time(ms) peri-stimulus')
ylabel('signed ERP')
title('Anterior Temporal Lobe')
% 
figure,
plot(taxis(tmin:tmax),D1(2,tmin:tmax,1),'r','LineWidth',3), hold on, plot(taxis(tmin:tmax),D1(2,tmin:tmax,2),'LineWidth',3)
% plot(D(2,500:950,1),'r'), hold on, plot(D(2,500:950,2))
xlabel('time(ms) peri-stimulus')
ylabel('signed ERP')
title('Inferior Frontal Gyrus')
% 
figure,
plot(taxis(tmin:tmax),D1(3,tmin:tmax,1),'r','LineWidth',3), hold on, plot(taxis(tmin:tmax),D1(3,tmin:tmax,2),'LineWidth',3)
% plot(D(3,500:950,1),'r'), hold on, plot(D(3,500:950,2))
xlabel('time(ms) peri-stimulus')
ylabel('signed ERP')
title('Middle Temporal Gyrus')

figure,
% plot(D(4,500:750,1),'r'), hold on, plot(D(4,500:750,2))
plot(taxis(tmin:tmax),D1(4,tmin:tmax,1),'r','LineWidth',3), hold on, plot(taxis(tmin:tmax),D1(4,tmin:tmax,2),'LineWidth',3)
xlabel('time(ms) peri-stimulus')
ylabel('signed ERP')
title('Angular Gyrus')
% 
figure,
plot(taxis(tmin:tmax),D1(5,tmin:tmax,1),'r','LineWidth',3), hold on, plot(taxis(tmin:tmax),D1(5,tmin:tmax,2),'LineWidth',3)
% plot(D(5,500:950,1),'r'), hold on, plot(D(5,500:950,2))
xlabel('time(ms) peri-stimulus')
ylabel('signed ERP')
title('Visual Word Form Area')