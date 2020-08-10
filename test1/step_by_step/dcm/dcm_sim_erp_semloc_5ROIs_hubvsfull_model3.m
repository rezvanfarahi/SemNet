clear
clc
close all
addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/')
% addpath('\\cbsu\data\imaging\local\software\spm_cbu_svn\releases\spm12_latest\')

% spm eeg
% spm('defaults','eeg')
% dbstop if error
% clear

Model = struct([]);



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

% LogEvd=[]; DCMname={};
this_model_dir='/imaging/rf02/TypLexMEG/dcm/5ROIs_hubvsfull/filtered_250ms_28models_5ROIs/semloc/maxCTF_ERP_dtr0/simulated/model3/';%'./16models_200ms/semloc';
out_dir='/imaging/rf02/TypLexMEG/dcm/5ROIs_hubvsfull/filtered_250ms_28models_5ROIs/semloc/maxCTF_ERP_dtr0/simulated/inverted/model3_stdSNR/';
if ~exist(out_dir,'dir')
    mkdir(out_dir)
end
parfor s = 1:length(dosubs) %parfor
    s
    sub    = dosubs(s);
    for n=[2,4,39:40,6,8,41:42,10,12,43:44,14,16,45:46,18:2:36,37:38]
        n
        this_filename=sprintf([out_dir,'DCM_simfin_erpf_SemLoc_5ROIs_maxCTF_sub%d_mod%d.mat'],sub,n);
        if ~exist(this_filename,'file')
            DCM1=load(sprintf([this_model_dir,'DCM_sim_erpf_SemLoc_5ROIs_maxCTF_sub%d_mod%d.mat'],sub,n));
            DCM=DCM1.DCM;
            DCM.name = sprintf([out_dir,'DCM_simfin_erpf_SemLoc_5ROIs_maxCTF_sub%d_mod%d.mat'],sub,n);
            
            %         DCM   = spm_dcm_erp(DCM);
            dcmcheck=0;
            while dcmcheck==0
                try
                    DCM   = spm_dcm_erp(DCM);
                    dcmcheck=1;
                end
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
