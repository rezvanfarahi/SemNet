clear
clc
close all
addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/')
% addpath('\\cbsu\data\imaging\local\software\spm_cbu_svn\releases\spm12_latest\')

spm eeg
spm('defaults','eeg')
cd('/imaging/rf02/TypLexMEG/dcm/latest/')

[FileName,PathName, filti] = uigetfile('*.mat','Select the DCM files', 'MultiSelect','on');
for ii=1:length(FileName)
    FileName{ii}
    load([PathName,FileName{ii}])%'/imaging/rf02/TypLexMEG/dcm/latest/150ms_7models/typlex/'
% figure,
spm_dcm_erp_results(DCM,'ERPs (mode)');%'ERPs (sources)''ERPs (mode)'
pause
% spm_dcm_erp_results(DCM,'Data')
end


