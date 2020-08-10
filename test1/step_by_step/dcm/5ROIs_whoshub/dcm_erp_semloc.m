clear
clc
close all
addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/')
% addpath('\\cbsu\data\imaging\local\software\spm_cbu_svn\releases\spm12_latest\')

spm eeg
spm('defaults','eeg')
% dbstop if error
% clear

Model = struct([]);

%% Create models in advance

% Columns = FROM, rows = TO

% For intrinsic (A) connections, {1} = forward; {2} = backward; {3} = lateral
Model(1).A{1} = [0 0 0 0 1
    0 0 0 0 1
    0 0 0 0 1
    0 0 0 0 1
    0 0 0 0 0];
Model(1).A{2} = [0 0 0 0 0 
    0 0 0 0 0 
    0 0 0 0 0 
    0 0 0 0 0 
    1 1 1 1 0];
Model(1).A{3} = [0 0 0 0 0 
    0 0 0 0 0 
    0 0 0 0 0 
    0 0 0 0 0 
    0 0 0 0 0];
%% ATL
% For modulatory (B) connections, {n} = nth effect   , bidirectional, WFA and ATL, input to WFA
Model(1).B{1} = [0 0 0 0 1 
    0 0 0 0 0 
    0 0 0 0 0 
    0 0 0 0 0 
    1 0 0 0 0];

Model(1).C =    [0
    0
    0
    0
    1];
% For modulatory (B) connections, {n} = nth effect   , WFA to ATL, input to WFA


Model(end+1) = Model(1);
Model(end).B{1} = [0 0 0 0 1 
    0 0 0 0 0 
    0 0 0 0 0 
    0 0 0 0 0 
    0 0 0 0 0];
Model(end).C =    [0
    0
    0
    0
    1];
% For modulatory (B) connections, {n} = nth effect   , ATL to WFA, input to WFA


Model(end+1) = Model(1);
Model(end).B{1} = [0 0 0 0 0 
    0 0 0 0 0 
    0 0 0 0 0 
    0 0 0 0 0 
    1 0 0 0 0];
Model(end).C =    [0
    0
    0
    0
    1];
%% SMG
% For modulatory (B) connections, {n} = nth effect   , bidirectional, WFA and SMG, input to WFA
Model(end+1) = Model(1);
Model(end).B{1} = [0 0 0 0 0 
    0 0 0 0 1 
    0 0 0 0 0 
    0 0 0 0 0 
    0 1 0 0 0];

Model(end).C =    [0
    0
    0
    0
    1];
% For modulatory (B) connections, {n} = nth effect   , WFA to SMG, input to WFA


Model(end+1) = Model(1);
Model(end).B{1} = [0 0 0 0 0 
    0 0 0 0 1 
    0 0 0 0 0 
    0 0 0 0 0 
    0 0 0 0 0];
Model(end).C =    [0
    0
    0
    0
    1];
% For modulatory (B) connections, {n} = nth effect   , SMG to WFA, input to WFA


Model(end+1) = Model(1);
Model(end).B{1} = [0 0 0 0 0 
    0 0 0 0 0 
    0 0 0 0 0 
    0 0 0 0 0 
    0 1 0 0 0];
Model(end).C =    [0
    0
    0
    0
    1];
%% IFG
% For modulatory (B) connections, {n} = nth effect   , bidirectional, WFA and IFG, input to WFA
Model(end+1) = Model(1);
Model(end).B{1} = [0 0 0 0 0 
    0 0 0 0 0 
    0 0 0 0 1 
    0 0 0 0 0 
    0 0 1 0 0];

Model(end).C =    [0
    0
    0
    0
    1];
% For modulatory (B) connections, {n} = nth effect   , WFA to IFG, input to WFA


Model(end+1) = Model(1);
Model(end).B{1} = [0 0 0 0 0 
    0 0 0 0 0 
    0 0 0 0 1 
    0 0 0 0 0 
    0 0 0 0 0];
Model(end).C =    [0
    0
    0
    0
    1];
% For modulatory (B) connections, {n} = nth effect   , IFG to WFA, input to WFA


Model(end+1) = Model(1);
Model(end).B{1} = [0 0 0 0 0 
    0 0 0 0 0 
    0 0 0 0 0 
    0 0 0 0 0 
    0 0 1 0 0];
Model(end).C =    [0
    0
    0
    0
    1];
%% MTG
% For modulatory (B) connections, {n} = nth effect   , bidirectional, WFA and MTG, input to WFA
Model(end+1) = Model(1);
Model(end).B{1} = [0 0 0 0 0 
    0 0 0 0 0 
    0 0 0 0 0 
    0 0 0 0 1 
    0 0 0 1 0];

Model(end).C =    [0
    0
    0
    0
    1];
% For modulatory (B) connections, {n} = nth effect   , WFA to MTG, input to WFA


Model(end+1) = Model(1);
Model(end).B{1} = [0 0 0 0 0 
    0 0 0 0 0 
    0 0 0 0 0 
    0 0 0 0 1 
    0 0 0 0 0];
Model(end).C =    [0
    0
    0
    0
    1];
% For modulatory (B) connections, {n} = nth effect   , ATL to MTG, input to WFA


Model(end+1) = Model(1);
Model(end).B{1} = [0 0 0 0 0 
    0 0 0 0 0 
    0 0 0 0 0 
    0 0 0 0 0 
    0 0 0 1 0];
Model(end).C =    [0
    0
    0
    0
    1];
%% ALL
% For modulatory (B) connections, {n} = nth effect   , bidirectional, WFA and ALL, input to WFA
Model(end+1) = Model(1);
Model(end).B{1} = [0 0 0 0 1 
    0 0 0 0 1 
    0 0 0 0 1 
    0 0 0 0 1 
    1 1 1 1 0];

Model(end).C =    [0
    0
    0
    0
    1];
% For modulatory (B) connections, {n} = nth effect   , WFA to the rest, input to WFA


Model(end+1) = Model(1);
Model(end).B{1} = [0 0 0 0 1 
    0 0 0 0 1 
    0 0 0 0 1 
    0 0 0 0 1 
    0 0 0 0 0];
Model(end).C =    [0
    0
    0
    0
    1];
% For modulatory (B) connections, {n} = nth effect   , ALL to WFA, input to WFA


Model(end+1) = Model(1);
Model(end).B{1} = [0 0 0 0 0 
    0 0 0 0 0 
    0 0 0 0 0 
    0 0 0 0 0 
    1 1 1 1 0];
Model(end).C =    [0
    0
    0
    0
    1];
% etc etc

% Null (baseline) model 1 - no modulation by Fac (input to both)

Model(end+1) = Model(1);
Model(end).B{1} = zeros(5);
% Model(end).C =    [1
%                  0];


%% Fit models

cwd = '/home/rf02/rezvan/test1/dcm/5ROIs_whoshub'

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
% Set up DCM structure

DCMbase = [];
DCMbase.xY.modality = 'LFP';
DCMbase.options.analysis = 'ERP';
DCMbase.options.model    = 'ERP';
DCMbase.options.spatial =  'LFP';
DCMbase.options.trials  = [1 2];  % [3 4] for CROI ([1 2] ROI)
DCMbase.options.Tdcm(1) = 0;      % start of peri-stimulus time to be modelled
DCMbase.options.Tdcm(2) = 200;   % end of peri-stimulus time to be modelled
DCMbase.options.Nmodes  = 8;      % nr of modes for data selection
DCMbase.options.h       = 0;      % nr of DCT components for detrending [increase to 4????]
DCMbase.options.onset   = 60;     % selection of onset (prior mean = 60ms for onset at 0)
DCMbase.options.D       = 1;      % downsampling
DCMbase.options.han     = 0;      % Hanning

DCMbase.options.lock     = 0;      % (if want modulations to be in same direction for all connections)
DCMbase.options.location = 0;      % (Not relevant; only for ECD)
DCMbase.options.symmetry = 0;      % (Not relevant; only for ECD)

DCMbase.options.nograph  = 0;

% Region names (Location, Lpos, irrelevant for LFP)
DCMbase.Sname = {'lATL','lSMG','lIFG','lMTG','lWFA'};
DCMbase.Lpos = [];

DCMbase.xU.X = [1; -1]; % W-NW
DCMbase.xU.name = {'C-A'};


% LogEvd=[]; DCMname={};
this_model_dir='./16models_200ms/semloc';
    if ~exist(this_model_dir,'dir')
        mkdir(this_model_dir)
    end
parfor s = 1:length(dosubs) %parfor
    
    sub    = dosubs(s);
    
    DCMsub = DCMbase;
    this_dir='./5ROIs_whoshub/';
    fname_a='_SemLoc_Evoked_5ROIs_forDCM';
    fname_in=[list_all{s},fname_a,'_cnvrt'];
    %     spm_eeg_load(fname_in)
    S=[]; S.D = sprintf(fname_in)
    
    DCMsub.xY.Dfile = S.D;
    
    DCMsub = spm_dcm_erp_data(DCMsub,DCMsub.options.h);  % Either convert your LFP data into an SPM object or replace this function with something that returns DCM structure with relevant field
    
    DCMsub = spm_dcm_erp_dipfit(DCMsub);  % Either convert your LFP data into an SPM object or replace this function with something that returns DCM structure with relevant field
    
    for n=1:numel(Model)
        s
        n
        DCM      = DCMsub;
        DCM.name = sprintf('./16models_200ms/semloc/DCM_erp_SemLoc_5ROIs_sub%d_mod%d.mat',sub,n);
        
        DCM.A = Model(n).A;
        DCM.B = Model(n).B;
        DCM.C = Model(n).C;
%         DCM   = spm_dcm_erp(DCM);
        dcmcheck=0;
        while dcmcheck==0
            try
                DCM   = spm_dcm_erp(DCM);
                dcmcheck=1;
            end
        end
        
        
        %         LogEvd(s,n) = DCM.F;
        %         DCMname{s,n} = DCM.name;
    end
end


%% Review results (have not cleaned)

% [dummy,maxi] = max(LogEvd);
%
% fprintf('Model %d is best',maxi);
%
% load(sprintf('DCM_erp_%d_%s',maxi,D.fname));
% % spm_dcm_erp_results_rik(DCM,'Coupling (C)');
% spm_dcm_erp_results_rik(DCM,'Coupling (A)');
% % spm_dcm_erp_results_rik(DCM,'Coupling (B)');
% spm_dcm_erp_results_rik(DCM,'Response');
%
% % spm_dcm_erp_results_rik(DCM,'ERPs (mode)');
% % spm_dcm_erp_results_rik(DCM,'ERPs (sources)');
%
%
% %    S=[]; S.D = sprintf('rtf_cbSPM_TFR_format_s%02d_retrieval_cr.mat',sub);
%     S=[]; S.D = sprintf('cbSPM_TFR_format_s%02d_retrieval_cr.mat',sub)
%
%     DCMsub.xY.Dfile = S.D;
%
%     DCMsub = spm_dcm_erp_dipfit(DCMsub, 1);
%
%     if ~isfield(DCMsub.xY,'source')  % (always true here)
%         DCMsub  = spm_dcm_ind_data(DCMsub);
%     end
%     DCMsub.options.gety = 0;
%     DCMsub.options.nograph  = 1;
%
%     for m=1:numel(model)
%
%         DCM      = DCMsub;
%         DCM.name = sprintf('DCM_ind_mod%d_sub%d_%s_%d-%d',m,sub,datafit,Twin(1),Twin(2));
%
%         DCM.A = model(m).A;
%         DCM.B = model(m).B;
%         DCM.C = model(m).C;
%
%         DCM   = spm_dcm_ind(DCM);   % saves automatically
%
%         LogEvd(s,m) = DCM.F;
%         DCMname{s,m} = DCM.name;
%     end
%
% end
%
