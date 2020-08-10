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
addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/')
% spm eeg
% spm('defaults','EEG');


% Some details about the data
%--------------------------------------------------------------------------
Nchannels = 2;
Nsamples  = 1201;
Ntrials   = 2;
TimeOnset = -0.5; % in sec
Fsample = 1000;

chlabels = {
            'lWFA' 
            'lATL'
        
            };

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
% define the output file name
for cnt=1:length(list_all)
%--------------------------------------------------------------------------
cnt

% create data array 
%--------------------------------------------------------------------------
fname_a='_SemLoc_Evoked_ATL_WFA_forDCM';
fname_in=[list_all{cnt},fname_a,'.mat'];
fname_out = [list_all{cnt},fname_a,'_cnvrt.mat'];
load(fname_in)
data = conabs_mat;%randn([Nchannels, Nsamples, Ntrials]);

% create the time axis (should be the same for all trials)
%--------------------------------------------------------------------------
timeaxis = [0:(Nsamples-1)]./Fsample + TimeOnset;
figure, subplot(2,1,1); plot(squeeze(data(1,:,:)))
subplot(2,1,2); plot(squeeze(data(2,:,:)))
end
