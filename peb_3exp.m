
clear

Models = [47, 2,4,39:40,6,8,41:42,10,12,43:44,14,16,45:46,18:2:36,37:38];  % 47 is full model

load('/imaging/rf02/Semnet/semnet4semloc/dcm/families_hierarchical_1st.mat')

fnames{1} = 'Full';
for f=1:length(family.names)
    fnames{f+1} = family.names{f};
end

%% Exp 1 (one task)
Subjs = [1:17]; Ns1 = length(Subjs) % Rezvan says ignore subject 1

% Get DCM.mat files
GCM1={};
for s=1:Ns1
    for m=1:length(Models) 
        tmp = load(sprintf('/imaging/rf02/Semnet/semnet4semloc/dcm/forPEB/250ms/SemLoc_SD/DCM_erpf_SemLoc_5ROIs_avg_maxCTF_sub%d_mod%d.mat',Subjs(s),Models(m)));
        GCM1{s,m} = tmp.DCM;
    end
end

% PEB
X = ones(Ns1,1);
PEB1 = spm_dcm_peb(GCM1, X);
spm_dcm_peb_review(PEB1,GCM1)

% BMR of full model
[BMA1,BMR1] = spm_dcm_peb_bmc(PEB1(1));
spm_dcm_peb_review(BMA1,GCM1)

% BMC across all 29 models
% Just create a 1D cell array that spm_dcm_peb_bmc needs
Mods = {}; for m = 1:length(Models); Mods{m} = GCM1{1,m}; end
[BMA1,BMR1] = spm_dcm_peb_bmc(PEB1(1),Mods);
spm_dcm_peb_review(BMA1,GCM1)

% BMC acros families
families = [1 family.partition+1]; % Need to add full model as new family!? (or see alternative below)
[BMA_fam1,fam1] = spm_dcm_peb_bmc_fam(BMA1, BMR1, families, 'WINNING');
[winP,winM] = max(fam1.family.post);
fprintf('Winning Family = %s (p=%4.3f)\n',fnames{winM},winP)

%This version excludes full model, though not sure my hacks to spm_dcm_peb_fam_rik totally valid
%families = [0 family.partition];
%[BMA_fam1_,fam] = spm_dcm_peb_bmc_fam_rik(BMA1, BMR1, families, 'WINNING');
%[winP,winM] = max(fam1.family.post);
%fprintf('Winning Family = %s (p=%4.3f)\n',family.names{winM},winP)


%% Exp2 (two tasks, within subject)
Subjs = [2:19]; Ns2 = length(Subjs) % Rezvan says ignore subject 1
GCM2a={};GCM2b={};
for s=1:Ns2
    for m=1:length(Models) 
        tmp = load(sprintf('/imaging/rf02/Semnet/semnet4semloc/dcm/forPEB/250ms/SemNet_LD/DCM_erpf_ConEmot_LD_oldreg_5ROIs_avg_maxCTF_sub%d_mod%d.mat',Subjs(s),Models(m)));
        GCM2a{s,m} = tmp.DCM;
        tmp = load(sprintf('/imaging/rf02/Semnet/semnet4semloc/dcm/forPEB/250ms/SemNet_SD/DCM_erpf_ConEmot_SL_oldreg_5ROIs_avg_maxCTF_sub%d_mod%d.mat',Subjs(s),Models(m)));
        GCM2b{s,m} = tmp.DCM;
    end
end

X = kron([1 1; 1 -1]',ones(Ns2,1));

GCM2 = [GCM2a; GCM2b];
PEB2 = spm_dcm_peb(GCM2, X);
spm_dcm_peb_review(PEB2,GCM2)

% BMR of full model
[BMA2,BMR2] = spm_dcm_peb_bmc(PEB2(1));
spm_dcm_peb_review(BMA2,GCM2)

% BMC across all 29 models
% Just create a 1D cell array that spm_dcm_peb_bmc needs
Mods = {}; for m = 1:length(Models); Mods{m} = GCM2{1,m}; end
[BMA2,BMR2] = spm_dcm_peb_bmc(PEB2(1),Mods);
spm_dcm_peb_review(BMA2,GCM2)

% BMC acros families
families = [1 family.partition+1]; % Need to add full model as new family!? (or see alternative below)
[BMA_fam2,fam2] = spm_dcm_peb_bmc_fam(BMA2, BMR2, families, 'WINNING');
Nfam = length(unique(families(find(families>0))))
[winP,winM] = max(fam2.family.post,[],[1 2],'linear');
[winMi,winMj] = ind2sub([Nfam Nfam],winM);
fprintf('Winning Family, Com = %s, Dif = %s, (p=%4.3f)\n',fnames{winMi},fnames{winMj},winP)

%This version excludes full model, though not sure my hacks to spm_dcm_peb_fam_rik totally valid
%families = [0 family.partition];
%[BMA_fam2_,fam2] = spm_dcm_peb_bmc_fam_rik(BMA2, BMR2, families, 'WINNING');
%Nfam = length(unique(families(find(families>0))))
%[winP,winM] = max(fam2.family.post,[],[1 2],'linear');
%[winMi,winMj] = ind2sub([Nfam Nfam],winM);
%fprintf('Winning Family, Com = %s, Dif = %s, (p=%4.3f)\n',fnames{winMi+1},fnames{winMj+1},winP)


%% All 3 exps in one PEB (though no modelling of within vs between subject error)

X = blkdiag(ones(Ns1,1), ones(Ns2,1), ones(Ns2,1));
X = X*[1 1 1; 1 -0.5 -0.5; 0 1 -1]'; % I don't think BMC includes 3rd-factor
X(:,2:3) = detrend(X(:,2:3),0);

GCM = [GCM1; GCM2a; GCM2b];
PEB = spm_dcm_peb(GCM, X);
spm_dcm_peb_review(PEB,GCM)

% BMR of full model
[BMA,BMR] = spm_dcm_peb_bmc(PEB(1));
spm_dcm_peb_review(BMA,GCM)

% BMC across all 29 models
% Just create a 1D cell array that spm_dcm_peb_bmc needs
Mods = {}; for m = 1:length(Models); Mods{m} = GCM{1,m}; end
[BMA,BMR] = spm_dcm_peb_bmc(PEB(1),Mods);
spm_dcm_peb_review(BMA,GCM)

% BMC acros families
families = [1 family.partition+1]; % Need to add full model as new family!? (or see alternative below)
[BMA_fam,fam] = spm_dcm_peb_bmc_fam(BMA, BMR, families, 'WINNING');
Nfam = length(unique(families(find(families>0))))
[winP,winM] = max(fam.family.post,[],[1 2],'linear');
[winMi,winMj] = ind2sub([Nfam Nfam],winM);
fprintf('Winning Family, Com = %s, Dif = %s, (p=%4.3f)\n',fnames{winMi},fnames{winMj},winP)

%This version excludes full model, though not sure my hacks to spm_dcm_peb_fam_rik totally valid
% families = [0 family.partition];
% [BMA_fam_,fam] = spm_dcm_peb_bmc_fam_rik(BMA, BMR, families, 'WINNING');
% Nfam = length(unique(families(find(families>0))))
% [winP,winM] = max(fam.family.post,[],[1 2],'linear');
% [winMi,winMj] = ind2sub([Nfam Nfam],winM);
% fprintf('Winning Family, Com = %s, Dif = %s, (p=%4.3f)\n',fnames{winMi+1},fnames{winMj+1},winP)


%% 3-level PEB (just for interest)

% Redo Exp2 PEB so that only commonalities (otherwise twice as many parameters as PEB1)
X = kron([1 1]',ones(Ns2,1));
GCM2 = [GCM2a; GCM2b];
PEB2 = spm_dcm_peb(GCM2, X);

PEBs = {PEB1(1); PEB2(1)} %{PEB1; PEB2};
X = [1 1; 1 -1]'
PEB3 = spm_dcm_peb(PEBs,X);

dF = [PEB.F] - PEB3.F;
figure,hist(dF)
onet(dF')



