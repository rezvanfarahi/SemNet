clear
clc
close all
addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/')
% addpath('\\cbsu\data\imaging\local\software\spm_cbu_svn\releases\spm12_latest\')

spm eeg
spm('defaults','eeg')

% clear

Model = struct([]);

%% Create models in advance

% Columns = FROM, rows = TO
% For intrinsic (A) connections, {1} = forward; {2} = backward; {3} = lateral
Model(1).A{1} = [0 0
    1 0];
Model(1).A{2} = [0 1
    0 0];
Model(1).A{3} = [0 0
    0 0];
% 1 For modulatory (B) connections, {n} = nth effect   , bidirectional, input to WFA
Model(1).B{1} = [0 1
    1 0];

Model(1).C =    [1
    0];
% 2 For modulatory (B) connections, {n} = nth effect   , bidirectional, input to ATL


Model(end+1) = Model(1);
Model(end).B{1} = [0 1
    1 0];
Model(end).C =    [0
    1];
% 3 modulatory (B) connections, {n} = nth effect , WFA to ATL, input to WFA

Model(end+1) = Model(1);
Model(end).B{1} = [0 0
    1 0];
Model(end).C =    [1
    0];

% 4 modulatory (B) connections, {n} = nth effect , WFA to ATL, input to ATL

Model(end+1) = Model(1);
Model(end).B{1} = [0 0
    1 0];
Model(end).C =    [0
    1];

% 5 modulatory (B) connections, {n} = nth effect , ATL to WFA, input to WFA

Model(end+1) = Model(1);
Model(end).B{1} = [0 1
    0 0];
Model(end).C =    [1
    0];
% 6 modulatory (B) connections, {n} = nth effect , ATL to WFA, input to ATL

Model(end+1) = Model(1);
Model(end).B{1} = [0 1
    0 0];
Model(end).C =    [0
    1];
% etc etc

% 7 Null (baseline) model 1 - no modulation by Fac (input to both)

Model(end+1) = Model(1);
Model(end).B{1} = zeros(2);
Model(end).C =    [1
    0];
% % 8 Null (baseline) model 1 - no modulation by Fac (input to ATL)
%
Model(end+1) = Model(1);
Model(end).B{1} = zeros(2);
Model(end).C =    [0
    1];
%% with self connections, modulation of WFA self-connection
% 9 For modulatory (B) connections, {n} = nth effect   , bidirectional, input to WFA
Model(end+1) = Model(1);

Model(end).B{1} = [1 1
    1 0];

Model(end).C =    [1
    0];

% 10 For modulatory (B) connections, {n} = nth effect   , bidirectional, input to ATL

%
Model(end+1) = Model(1);
Model(end).B{1} = [1 1
    1 0];
Model(end).C =    [0
    1];
% 11 modulatory (B) connections, {n} = nth effect , WFA to ATL, input to WFA

Model(end+1) = Model(1);
Model(end).B{1} = [1 0
    1 0];
Model(end).C =    [1
    0];

% 12 modulatory (B) connections, {n} = nth effect , WFA to ATL, input to ATL

Model(end+1) = Model(1);
Model(end).B{1} = [1 0
    1 0];
Model(end).C =    [0
    1];

% 13 modulatory (B) connections, {n} = nth effect , ATL to WFA, input to WFA

Model(end+1) = Model(1);
Model(end).B{1} = [1 1
    0 0];
Model(end).C =    [1
    0];
% 14 modulatory (B) connections, {n} = nth effect , ATL to WFA, input to ATL

Model(end+1) = Model(1);
Model(end).B{1} = [1 1
    0 0];
Model(end).C =    [0
    1];

% etc etc

% 15 Null (baseline) model 1 - no modulation by Fac (input to WFA)

Model(end+1) = Model(1);
Model(end).B{1} = [1 0
    0 0];
Model(end).C =    [1
    0];

% 16 Null (baseline) model 1 - no modulation by Fac (input to ATL)

Model(end+1) = Model(1);
Model(end).B{1} = [1 0
    0 0];
Model(end).C =    [0
    1];

%% with self connections- modulation of ATL self connections

% 17 For modulatory (B) connections, {n} = nth effect   , bidirectional, input to WFA
Model(end+1) = Model(1);
Model(end).B{1} = [0 1
    1 1];
Model(end).C =    [1
    0];
% 18 For modulatory (B) connections, {n} = nth effect   , bidirectional, input to ATL

%
Model(end+1) = Model(1);
Model(end).B{1} = [0 1
    1 1];
Model(end).C =    [0
    1];
% 19 modulatory (B) connections, {n} = nth effect , WFA to ATL, input to WFA

Model(end+1) = Model(1);
Model(end).B{1} = [0 0
    1 1];
Model(end).C =    [1
    0];

% 20 modulatory (B) connections, {n} = nth effect , WFA to ATL, input to ATL

Model(end+1) = Model(1);
Model(end).B{1} = [0 0
    1 1];
Model(end).C =    [0
    1];

% 21 modulatory (B) connections, {n} = nth effect , ATL to WFA, input to WFA

Model(end+1) = Model(1);
Model(end).B{1} = [0 1
    0 1];
Model(end).C =    [1
    0];
% 22 modulatory (B) connections, {n} = nth effect , ATL to WFA, input to ATL

Model(end+1) = Model(1);
Model(end).B{1} = [0 1
    0 1];
Model(end).C =    [0
    1];

% etc etc

% 23 Null (baseline) model 1 - no modulation by Fac (input to WFA)

Model(end+1) = Model(1);
Model(end).B{1} = [0 0
    0 1];
Model(end).C =    [1
    0];

% 24 Null (baseline) model 1 - no modulation by Fac (input to ATL)

Model(end+1) = Model(1);
Model(end).B{1} = [0 0
    0 1];
Model(end).C =    [0
    1];

%% with self connections- modulation of both ATL and WFA self connections

% 25 For modulatory (B) connections, {n} = nth effect   , bidirectional, input to WFA
Model(end+1) = Model(1);
Model(end).B{1} = [1 1
    1 1];
Model(end).C =    [1
    0];
% 26 For modulatory (B) connections, {n} = nth effect   , bidirectional, input to ATL

%
Model(end+1) = Model(1);
Model(end).B{1} = [1 1
    1 1];
Model(end).C =    [0
    1];
% 27 modulatory (B) connections, {n} = nth effect , WFA to ATL, input to WFA

Model(end+1) = Model(1);
Model(end).B{1} = [1 0
    1 1];
Model(end).C =    [1
    0];

% 28 modulatory (B) connections, {n} = nth effect , WFA to ATL, input to ATL

Model(end+1) = Model(1);
Model(end).B{1} = [1 0
    1 1];
Model(end).C =    [0
    1];

% 29 modulatory (B) connections, {n} = nth effect , ATL to WFA, input to WFA

Model(end+1) = Model(1);
Model(end).B{1} = [1 1
    0 1];
Model(end).C =    [1
    0];
% 30 modulatory (B) connections, {n} = nth effect , ATL to WFA, input to ATL

Model(end+1) = Model(1);
Model(end).B{1} = [1 1
    0 1];
Model(end).C =    [0
    1];

% etc etc

% 31 Null (baseline) model 1 - no modulation by Fac (input to WFA)

Model(end+1) = Model(1);
Model(end).B{1} = [1 0
    0 1];
Model(end).C =    [1
    0];

% 32 Null (baseline) model 1 - no modulation by Fac (input to ATL)

Model(end+1) = Model(1);
Model(end).B{1} = [1 0
    0 1];
Model(end).C =    [0
    1];

%% Fit models

cwd = '/home/rf02/rezvan/test1/dcm'

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
DCMbase.options.Tdcm(2) = 150;   % end of peri-stimulus time to be modelled
DCMbase.options.Nmodes  = 8;      % nr of modes for data selection
DCMbase.options.h       = 0;      % nr of DCT components for detrending [increase to 4????]
DCMbase.options.onset   = 60;     % selection of onset (prior mean = 60ms for onset at 0)
DCMbase.options.D       = 1;      % downsampling
DCMbase.options.han     = 0;      % Hanning

DCMbase.options.lock     = 0;      % (if want modulations to be in same direction for all connections)
DCMbase.options.location = 0;      % (Not relevant; only for ECD)
DCMbase.options.symmetry = 0;      % (Not relevant; only for ECD)

DCMbase.options.nograph  = 0;
% DCMbase.options.Nmax    = 128;



% Region names (Location, Lpos, irrelevant for LFP)
DCMbase.Sname = {'lWFA','lATL'};
DCMbase.Lpos = [];

DCMbase.xU.X = [1; -1]; % W-NW
DCMbase.xU.name = {'W-NW'};
% main_dir='/imaging/rf02/TypLexMEG/dcm/latest/';
this_dir='/imaging/rf02/TypLexMEG/dcm/latest/150ms_32models_2ROIs/typlex';%'/home/rf02/rezvan/test1/step_by_step/dcm/latest/200ms_7models/typlex';
    if ~exist(this_dir,'dir')
        mkdir(this_dir)
    end
LogEvd=[]; DCMname={};

parfor s = 1:length(dosubs) %parfor
    
    sub    = dosubs(s);
    
    DCMsub = DCMbase;
    fname_a='_TypLex_Evoked_2ROIs_forDCM';%'_TypLex_Evoked_ATL_WFA_forDCM';
    path_in='/imaging/rf02/TypLexMEG/dcm/2ROIs_whoshub/';
    fname_in=[path_in,list_all{sub}(3:end),fname_a,'_cnvrt'];
%     load(fname_in)
    S=[]; S.D = sprintf(fname_in)
    
    DCMsub.xY.Dfile = S.D;
%     DCMsub.xY.Ic=[5,4];
    
    DCMsub = spm_dcm_erp_data(DCMsub,DCMsub.options.h);  % Either convert your LFP data into an SPM object or replace this function with something that returns DCM structure with relevant field
    
    DCMsub = spm_dcm_erp_dipfit(DCMsub);  % Either convert your LFP data into an SPM object or replace this function with something that returns DCM structure with relevant field
    
    for n=[1,8,16,17,25]%1:numel(Model)
        sub
        n
%         nnew=n+20;
        DCM      = DCMsub;
        DCM.name = sprintf([this_dir,'/DCM_erp_TypLex_ATL_WFA_new_sub%d_mod%d.mat'],sub,n);
        
        DCM.A = Model(n).A;
        DCM.B = Model(n).B;
        DCM.C = Model(n).C;
        
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
