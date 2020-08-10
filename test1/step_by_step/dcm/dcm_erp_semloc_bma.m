clear 
clc
addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r6470')%addpath('/imaging/rf02/TypLexMEG/dcm/spm_12_latest/')
spm eeg

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
LogEvd=[]; DCMname={};
all_DCMs={};
subj=struct();
GCM={};
dcm_path='/imaging/rf02/TypLexMEG/dcm/5ROIs_hubvsfull/filtered_450ms_28models_5ROIs/semloc/maxCTF_ERP_dtr0/';%'/imaging/rf02/TypLexMEG/dcm/5ROIs_hubvsfull/250ms_38models_5ROIs/semloc0/';
for s = 1:length(dosubs) %parfor
sess=struct();
sub    = dosubs(s);
ncnt=0;
for n=[14,16,45,46]%[2,4,39,40]%,6,8,41:42,10,12,43:44,14,16,45:46,18:2:36,37:38]%[2:2:36,37:38]%1:38%numel(Model)
ncnt=ncnt+1;
sub
n
load(sprintf([dcm_path,'DCM_erpf_SemLoc_5ROIs_avg_maxCTF_sub%d_mod%d.mat'],sub,n));
GCM{s,ncnt}=DCM;
end
end
%% rfx
% [n,m] = size(GCM);
% for i = 1:n
% for j = 1:m
% if ~isfield(GCM{i,j}, 'Ep')
% error(['Could not average: subject %d model %d ' ...
% 'not estimated'], i, j);
% end
% subj(i).sess(1).model(j).Ep = GCM{i,j}.Ep;
% subj(i).sess(1).model(j).Cp = GCM{i,j}.Cp;
% F(i,j) = GCM{i,j}.F;
% subj(i).sess(1).model(j).fname=GCM{i,j}.name;
% end
% end
% F    = sum(F,1);
% F    = F - repmat(max(F')',1,m);
% P    = exp(F);
% post = P./repmat(sum(P,2),1,m);
% indx = 1:m;
% bma        = spm_dcm_bma(post,indx,subj,10000);
%% ffx
[n,m] = size(GCM);
for i = 1:n
for j = 1:m
if ~isfield(GCM{i,j}, 'Ep')
error(['Could not average: subject %d model %d ' ...
'not estimated'], i, j);
end
subj(i).sess(1).model(j).Ep = GCM{i,j}.Ep;
subj(i).sess(1).model(j).Cp = GCM{i,j}.Cp;
F(i,j) = GCM{i,j}.F;
end
end
F    = sum(F,1);
F    = F - max(F);
P    = exp(F);
post = P/sum(P);
indx = 1:m;
bma        = spm_dcm_bma(GCM);%(post,indx,subj,10000);

a=bma.Ep.B{1,1}./(bma.Cp.B{1,1}/sqrt(17));
% b=a(~isnan(a));
% v=16;
% p=[];
% for i=1:length(b)
%     t=b(i);
% p(i)=(1-betainc(v/(v+t^2),v/2,0.5));
% end
thispar=zeros(5,5,17);
for s=1:17
%     for ii=1:5
%         for jj=1:5
thispar(:,:,s)=bma.SUB(s).Ep.B{1,1};%bma.mEps{1,s}.B{1,1};
%         end
%     end
end
p=[];
    for ii=1:5
        for jj=1:5
            [h,p(ii,jj)]=ttest(thispar(ii,jj,:));
        end
    end
thisparm=squeeze(mean(thispar,3));
thispars=squeeze(std(thispar,0,3));
a=thisparm./(thispars/sqrt(17));
% bma.mEp.B{1,1}./(bma.sEp.B{1,1}/sqrt(17))
% a=ans;
% b=a(~isnan(a));
% b=1-tcdf(b,16);

Fmat=[];
for ii=1:17
for jj=1:4
Fmat(ii,jj)=GCM{ii,jj}.F;
end
end
mean(Fmat,1);

parall=zeros(5,5,17);
for ii=1:17
parall(:,:,ii)=GCM{ii,4}.Ep.B{1,1};
end
parmean=mean(parall,3);


parallmean=zeros(5,5,17);
aa=randperm(17);
cc=0;
for ii=aa
    cc=cc+1;
parallmean(:,:,cc)=squeeze(mean(parall(:,:,aa(1:cc)),3));
end
cc=0;
for jj=aa
    cc=cc+1;
thisbpa= spm_dcm_bpa(GCM(aa(1:cc),4),1);
bpamean(:,:,cc)=thisbpa.Ep.B{1,1};
end