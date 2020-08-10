clear
clc
close all
addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/')
% addpath('\\cbsu\data\imaging\local\software\spm_cbu_svn\releases\spm12_latest\')

% spm eeg
% spm('defaults','eeg')
% 
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
% For modulatory (B) connections, {n} = nth effect   , bidirectional, input to WFA
Model(1).B{1} = [0 1
                1 0];

Model(1).C =    [1
                 0];
% For modulatory (B) connections, {n} = nth effect   , bidirectional, input to ATL

% 
Model(end+1) = Model(1);
Model(end).B{1} = [0 1
                   1 0];
Model(end).C =    [0
                   1];
% modulatory (B) connections, {n} = nth effect , WFA to ATL, input to WFA

Model(end+1) = Model(1);
Model(end).B{1} = [0 0
                 1 0];
Model(end).C =    [1
                 0];
             
% modulatory (B) connections, {n} = nth effect , WFA to ATL, input to ATL

Model(end+1) = Model(1);
Model(end).B{1} = [0 0
                 1 0];
Model(end).C =    [0
                 1];

% modulatory (B) connections, {n} = nth effect , ATL to WFA, input to WFA

Model(end+1) = Model(1);
Model(end).B{1} = [0 1
                 0 0];
Model(end).C =    [1
                 0];
% modulatory (B) connections, {n} = nth effect , ATL to WFA, input to ATL

Model(end+1) = Model(1);
Model(end).B{1} = [0 1
                 0 0];
Model(end).C =    [0
                 1];             
% etc etc

% Null (baseline) model 1 - no modulation by Fac (input to both)

Model(end+1) = Model(1);
Model(end).B{1} = zeros(2);
% Model(end).C =    [1
%                   1]; 


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



% Region names (Location, Lpos, irrelevant for LFP)
DCMbase.Sname = {'lWFA','lATL'};
DCMbase.Lpos = [];

DCMbase.xU.X = [1; -1]; % W-NW
DCMbase.xU.name = {'W-NW'};

LogEvd=[]; DCMname={};

for s = 1:length(dosubs) %parfor
    
    sub    = dosubs(s);
    
   
    for n=1:numel(Model)
        sub
        n

        load(sprintf('./latest/150ms_7models/typlex/DCM_erp_TypLex_ATL_WFA_sub%d_mod%d.mat',sub,n));
       
        
        LogEvd(s,n) = DCM.F;
%         DCMname{s,n} = DCM.name;
    end
end
[alpha,exp_r,xp,pxp,bor] = spm_BMS (LogEvd, [], 1, [], [], []);
% for n1=1:numel(Model)
%     subplot(3,3,n1), bar(LogEvd(1:17,n1), 'DisplayName', 'LogEvd', 'YDataSource', 'LogEvd(1:17,n1)'); figure(gcf)
% end

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
