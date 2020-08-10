clear
clc
% close all
addpath('/imaging/rf02/TypLexMEG/dcm/spm_12_latest/')%('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/')%('/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r6685/')%
% addpath('\\cbsu\data\imaging\local\software\spm_cbu_svn\releases\spm12_latest\')

spm eeg
% spm('defaults','eeg')
% 
% clear




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
which_model=13;
out_dir='/imaging/rf02/TypLexMEG/dcm/5ROIs_hubvsfull/filtered_250ms_28models_5ROIs/semloc/maxCTF_ERP_dtr0/simulated_thesis/model13/';%'/home/rf02/rezvan/test1/step_by_step/dcm/latest/200ms_7models/semloc';
if ~exist(out_dir,'dir')
    mkdir(out_dir)
end

% Region names (Location, Lpos, irrelevant for LFP)
% DCMbase.Sname = {'lWFA','lATL'};
% DCMbase.Lpos = [];
% 
% DCMbase.xU.X = [1; -1]; % W-NW
% DCMbase.xU.name = {'W-NW'};

LogEvd=[]; DCMname={};
all_DCMs={};
subj=struct();
GCM={};
% GCM_final={};

dcm_path='/imaging/rf02/TypLexMEG/dcm/5ROIs_hubvsfull/filtered_250ms_28models_5ROIs/semloc/maxCTF_ERP_dtr0/';%'/imaging/rf02/TypLexMEG/dcm/5ROIs_hubvsfull/250ms_38models_5ROIs/semloc0/';

for s = 1:length(dosubs) %parfor
    sess=struct();
    sub    = dosubs(s);    
   ncnt=0;
    for n=[2,4,39:40,6,8,41:42,10,12,43:44,14,16,45:46,18:2:36,37:38]%[2:2:36,37:38]%1:38%numel(Model)
        ncnt=ncnt+1;
        sub
        n
        load(sprintf([dcm_path,'DCM_erpf_SemLoc_5ROIs_avg_maxCTF_sub%d_mod%d.mat'],sub,n));
        GCM{s,ncnt}=DCM;
    end
end
% for mcnt=1:size(GCM,2)%no of models
%     [GCMpre,gen] = spm_dcm_simulate_latestVladimir(GCM, 'estimated', 1, mcnt);
%     GCM_final(:,mcnt)=GCMpre(:,mcnt);
% end
% GCM_final={};
% for subno=1:17
% [GCM_final(subno,:),gen] = spm_dcm_simulate_latestVladimir(GCM(subno,:), 'var', [], which_model);
% end
[GCM_final,gen] = spm_dcm_simulate_latestVladimir(GCM, 'estimated', [], which_model);

nsubs=size(GCM,1);
nmodls=size(GCM,2);
% for gcnt1=1:nsubs
%     for gcnt2=1:nmodls
%         GCM_final{gcnt1,gcnt2}.options.DATA=1;
%     end
% end


GCM_meannsqr=zeros(nsubs,nmodls,2);
for gcnt1=1:nsubs
    for gcnt2=1:nmodls
        for tcnt=1:2 %no of conditions
            meansqrt=(GCM{gcnt1, gcnt2}.xY.y{1,tcnt}-GCM_final{gcnt1, gcnt2}.xY.y{1,tcnt}).^2;
            meansqrt1=sqrt(mean(meansqrt,1));
            GCM_meannsqr(gcnt1,gcnt2,tcnt)=mean(meansqrt1);
        end
    end
end
GCM_meannsqrf=mean(GCM_meannsqr,3);
% save([out_dir,'GCM_meannsqrf_var',num2str(which_model),'.mat'],'GCM_meannsqrf')
figure, boxplot(GCM_meannsqrf(:,1))
mean(GCM_meannsqrf(:))
GCM_hub=mean(GCM_meannsqrf(:,1:16),2);
GCM_nohub=mean(GCM_meannsqrf(:,17:26),2);
GCM_nomodu=mean(GCM_meannsqrf(:,27:end),2);
GCM_hubvselse=[GCM_hub,GCM_nohub,GCM_nomodu];
mean(GCM_hub(:))
mean(GCM_nohub(:))
mean(GCM_nomodu(:))

GCM_atlhub=mean(GCM_meannsqrf(:,1:4),2);
GCM_ifghub=mean(GCM_meannsqrf(:,5:8),2);
GCM_mtghub=mean(GCM_meannsqrf(:,9:12),2);
GCM_aghub=mean(GCM_meannsqrf(:,13:16),2);
GCM_hubvselse=[GCM_atlhub,GCM_ifghub,GCM_mtghub,GCM_aghub];
mean(GCM_atlhub(:))
mean(GCM_ifghub(:))
mean(GCM_mtghub(:))
mean(GCM_aghub(:))

GCM_atlfam=mean(GCM_meannsqrf(:,1:4),1);
GCM_atlfam


subj=struct();
parfor s = 1:size(GCM_final,1) %parfor
    sess=struct();
    sub    = dosubs(s);
    
    ncnt=0;
    for n=[2,4,39:40,6,8,41:42,10,12,43:44,14,16,45:46,18:2:36,37:38]%[2:2:36,37:38]%1:38%numel(Model)
        ncnt=ncnt+1;
        sub
        n       
        DCM1=GCM_final{s,ncnt};%[GCM_final(subno,:),gen] = spm_dcm_simulate_latestVladimir(GCM(subno,:), 'var', [], which_model);

        DCM1.name = sprintf([out_dir,'/DCM_simfin_erpf_SemLoc_5ROIs_maxCTF_sub%d_mod%d.mat'],sub,n);
        %         DCM   = spm_dcm_erp(DCM);
        dcmcheck=0;
        while dcmcheck==0
            try
                DCM   = spm_dcm_erp(DCM1);
                dcmcheck=1;
            end
        end
%         save(DCM.name, 'DCM', spm_get_defaults('mat.format'));
%         load(sprintf([out_dir,'DCM_sim_erpf_SemLoc_5ROIs_maxCTF_sub%d_mod%d.mat'],sub,n));
%         this_model=struct();
%         this_model.fname=sprintf([out_dir,'DCM_sim_erpf_SemLoc_5ROIs_maxCTF_sub%d_mod%d.mat'],sub,n);%DCM.name;
%         this_model.F=DCM.F;
% 
%         this_model.Ep=DCM.Ep;
%         this_model.Cp=DCM.Cp;
%         sess.model(1,ncnt)=this_model;
%         LogEvd(s,ncnt) = DCM.F;
%         DCMname{s,ncnt} = DCM.name;
%         all_DCMs{s,ncnt}=DCM;
    end
%     subj(1,s).sess=sess;
end

% save([out_dir,'simulated_model1_space_allmodels.mat'],'subj')

% save(DCM.name, 'DCM', spm_get_defaults('mat.format'));
