clear
clc
close all
% addpath('/imaging/rf02/TypLexMEG/dcm/spm_12_latest/')
addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r6685/')%('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/')%
% addpath('\\cbsu\data\imaging\local\software\spm_cbu_svn\releases\spm12_latest\')

spm eeg
% spm('defaults','eeg')
% 
% clear

Model = struct([]);

%% Create models in advance

% % Columns = FROM, rows = TO
% % For intrinsic (A) connections, {1} = forward; {2} = backward; {3} = lateral
% Model(1).A{1} = [0 0
%                  1 0];
% Model(1).A{2} = [0 1
%                 0 0];
% Model(1).A{3} = [0 0
%                 0 0];
% % For modulatory (B) connections, {n} = nth effect   , bidirectional, input to WFA
% Model(1).B{1} = [0 1
%                 1 0];
% 
% Model(1).C =    [1
%                  0];
% % For modulatory (B) connections, {n} = nth effect   , bidirectional, input to ATL
% 
% % 
% Model(end+1) = Model(1);
% Model(end).B{1} = [0 1
%                    1 0];
% Model(end).C =    [0
%                    1];
% % modulatory (B) connections, {n} = nth effect , WFA to ATL, input to WFA
% 
% Model(end+1) = Model(1);
% Model(end).B{1} = [0 0
%                  1 0];
% Model(end).C =    [1
%                  0];
%              
% % modulatory (B) connections, {n} = nth effect , WFA to ATL, input to ATL
% 
% Model(end+1) = Model(1);
% Model(end).B{1} = [0 0
%                  1 0];
% Model(end).C =    [0
%                  1];
% 
% % modulatory (B) connections, {n} = nth effect , ATL to WFA, input to WFA
% 
% Model(end+1) = Model(1);
% Model(end).B{1} = [0 1
%                  0 0];
% Model(end).C =    [1
%                  0];
% % modulatory (B) connections, {n} = nth effect , ATL to WFA, input to ATL
% 
% Model(end+1) = Model(1);
% Model(end).B{1} = [0 1
%                  0 0];
% Model(end).C =    [0
%                  1];             
% % etc etc
% 
% % Null (baseline) model 1 - no modulation by Fac (input to both)
% 
% Model(end+1) = Model(1);
% Model(end).B{1} = zeros(2);
% % Model(end).C =    [1
% %                   1]; 


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
% DCMbase.Sname = {'lWFA','lATL'};
% DCMbase.Lpos = [];
% 
% DCMbase.xU.X = [1; -1]; % W-NW
% DCMbase.xU.name = {'W-NW'};

LogEvd=[]; DCMname={};
Amean=zeros(length(dosubs),28,5,5,2);
Bmean=zeros(length(dosubs),28,5,5);
all_DCMs={};
subj=struct();
% dcm_path='/imaging/rf02/TypLexMEG/dcm/5ROIs_hubvsfull/filtered_250ms_28models_5ROIs/semloc/maxCTF_ERP_dtr0/simulated/inverted/model1/';%'/imaging/rf02/TypLexMEG/dcm/5ROIs_hubvsfull/filtered_250ms_28models_5ROIs/semloc/maxCTF_aifg_ERP_dtr1/';%'/imaging/rf02/TypLexMEG/dcm/5ROIs_hubvsfull/filtered_250ms_28models_5ROIs/semloc/avg_allverts_dtr0_dip/dip/';%/simulated/inverted/';%'/imaging/rf02/TypLexMEG/dcm/5ROIs_hubvsfull/filtered_250ms_28models_5ROIs/semloc/';
% dcm_path='/imaging/rf02/TypLexMEG/dcm/5ROIs_hubvsfull/filtered_250ms_28models_5ROIs/semloc/maxCTF_ERP_dtr0/simulated_thesis/model15/';
dcm_path='/imaging/rf02/TypLexMEG/dcm/5ROIs_hubvsfull/filtered_250ms_28models_5ROIs/semloc/maxCTF_ERP_dtr0_fixedA0/';%'/imaging/rf02/TypLexMEG/dcm/5ROIs_hubvsfull/filtered_250ms_28models_5ROIs/semloc/maxCTF_ERP_dtr0/';

% dcm_path='/imaging/rf02/TypLexMEG/dcm/5ROIs_hubvsfull/filtered_250ms_28models_5ROIs/semloc/maxCTF_ERP_dtr0/simulated/inverted/model3/';


% /imaging/rf02/TypLexMEG/dcm/5ROIs_hubvsfull/filtered_250ms_28models_5ROIs/semloc/avg_50verts_dtr1_dip/ctf

for s = 1:length(dosubs) %parfor
    sess=struct();
    sub    = dosubs(s);
    
   ncnt=0;
    for n=[2,4,39:40,6,8,41:42,10,12,43:44,14,16,45:46,18:2:36,37:38]%[2:2:36,37:38]%1:38%numel(Model)
        ncnt=ncnt+1;
        sub
        n

        load(sprintf([dcm_path,'DCM_erpf_SemLoc_5ROIs_avg_maxCTF_sub%d_mod%d.mat'],sub,n));%'DCM_erpf_SemLoc_5ROIs_sub%d_mod%d.mat'
       
        this_model=struct();
this_model.fname=sprintf([dcm_path,'DCM_erpf_SemLoc_5ROIs_avg_maxCTF_sub%d_mod%d.mat'],sub,n);%DCM.name;DCM_erpf_SemLoc_5ROIs_mCTFavg_sub14_mod38
this_model.F=DCM.F;

this_model.Ep=DCM.Ep;
this_model.Cp=DCM.Cp;
sess.model(1,ncnt)=this_model;
        LogEvd(s,ncnt) = DCM.F;
        DCMname{s,ncnt} = DCM.name;
        all_DCMs{s,ncnt}=DCM;
        Bmean(s,ncnt,:,:)=DCM.Ep.B{1,1};
        Amean(s,ncnt,:,:,1)=DCM.Ep.A{1,1};
        Amean(s,ncnt,:,:,2)=DCM.Ep.A{1,2};
        
    end
    subj(1,s).sess=sess;
end
out_dir=dcm_path;%'/imaging/rf02/TypLexMEG/dcm/5ROIs_hubvsfull/filtered_250ms_28models_5ROIs/semloc/avg_allverts_dtr0_dip/dip/';%/simulated/inverted/';%'/home/rf02/rezvan/test1/step_by_step/dcm/latest/200ms_7models/semloc';
if ~exist(out_dir,'dir')
    mkdir(out_dir)
end
save([out_dir,'/modelspace_allmaxCTF_all_f250_fixedA0.mat'],'subj')
cd(dcm_path)
% [alpha,exp_r,xp,pxp,bor] = spm_BMS (LogEvd, [], 1, [], [], []);
% for n1=1:numel(Model)
%     subplot(3,3,n1), bar(LogEvd(1:17,n1), 'DisplayName', 'LogEvd', 'YDataSource', 'LogEvd(1:17,n1)'); figure(gcf)
% end
[post,exp_r,xp,pxp,bor] = spm_dcm_bmc(all_DCMs);
% Bayesian model comparison
% FORMAT [post,exp_r,xp,pxp,bor] = spm_dcm_bmc(DCM)
%
% DCM     - {subjects x models} cell array of DCMs
% ------------------------------------------------
%     DCM{i,j}.F  - free energy
%
% OUTPUTS
% -------
% post    - FFX posterior model probabilities p(m|y)
% exp_r   - RFX expectation of the posterior  p(m|y)
% xp      - RFX exceedance probabilities
% pxp     - RFX protected exceedance probabilities
% bor     - RFX Bayes Omnibus Risk (probability that model frequencies 
%           are equal)
%
% This routine computes fixed and random effects posterior probabilities
% over models. It also returns exceedance  probabilities and protected
% statistics.
% 
% See also: spm_dcm_bma.m and spm_BMS.m


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
