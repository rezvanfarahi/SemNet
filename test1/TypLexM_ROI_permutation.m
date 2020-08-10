clear
clc
addpath('/home/rf02/rezvan/test1/BCT/2015_01_25 BCT')
in_path='/imaging/rf02/TypLexMEG/';

list_all = {'meg10_0378/101209/',
    'meg10_0390/101214/',
    'meg11_0026/110223/',
    'meg11_0050/110307/',
    'meg11_0052/110307/',
    'meg11_0069/110315/',
    'meg11_0086/110322/',
    'meg11_0091/110328/',
    'meg11_0096/110404/',
    'meg11_0101/110411/',
    'meg11_0102/110411/',
    'meg11_0112/110505/',
    'meg11_0104/110412/',
    'meg11_0118/110509/',
    'meg11_0131/110519/',
    'meg11_0144/110602/',
    'meg11_0147/110603/',
    };

fname_a='pca_ppc_abstract_cwtmorlet_resmpled.mat';
fname_c='pca_ppc_concrete_cwtmorlet_resmpled.mat';
fname_b='pca_ppc_bothcon_cwtmorlet_resmpled.mat';

fname=[in_path,list_all{1},fname_a];
load(fname)
% fname_s='plv_subtract_cwtmorlet.mat';
n_subj=length(list_all);
node_perc=[0.05:0.05:0.75];

tr=zeros(n_subj,1);
trb=zeros(n_subj,1);

DEG_br=zeros(46,4,2,17);
load('/imaging/rf02/TypLexMEG/icaanalysis_results/jpg/time_resolved/ROI_connectivity/graph/labels_46.mat')
load('/imaging/rf02/TypLexMEG/icaanalysis_results/jpg/time_resolved/ROI_connectivity/graph/node_angles46.mat')
for cnt=1:n_subj

    load([in_path,list_all{cnt},'nodestrength_15_pca_ppc_concrete_resampled_2win.mat'])
    load([in_path,list_all{cnt},'nodestrength_15_pca_ppc_abstract_resampled_2win.mat'])
    deg_cr=strength_c;
    deg_ar=strength_a;

%     load([in_path,list_all{cnt},'degrees_proportional_15_pca_ppc_both_resampled.mat'])
    degb=cat(4,deg_cr,deg_ar);deg_br=mean(degb,4);
    DEG_br(:,:,:,cnt)=deg_cr-deg_ar;
end
load('/imaging/rf02/TypLexMEG/icaanalysis_results/jpg/time_resolved/ROI_connectivity/graph/ppc_degrees_15_resampled_2win_bothcon2_hubs_prob.mat');

percheckp=zeros(1,5000);
percheckn=zeros(1,5000);
for percnt=1:5000
    percnt
    deg_brt=zeros(size(DEG_br));
    rcoef=sign(randn(1,17));
    for cntt=1:17
        if rcoef(cntt)~=0
deg_brt(:,:,:,cntt)=rcoef(cntt)*DEG_br(:,:,:,cntt);%(:,:,40:end,:);
        else
            deg_brt(:,:,:,cntt)=sign(randn(1,1))*DEG_br(:,:,:,cntt);%(:,:,40:end,:);
        end
    end
%deg_brr=DEG_br(:,:,1:40,:);


deg_brtr=deg_brt;%(:,:,2:2:10,:);%zeros(size(deg_brt,1),size(deg_brt,2),10,size(deg_brt,4));

deg_brtzr=zeros(size(deg_brtr));


    for kk=1:size(deg_brtr,4)%subjs
        deg_brtzr(:,:,:,kk)=zscore(deg_brtr(:,:,:,kk));
%         windeg1=deg_brr(:,jj,:,kk);
%         deg_brrz(:,jj,:,kk)=(deg_brr(:,jj,:,kk)-mean(windeg1(:)))/std(windeg1(:));
    end
   

deg_brtzrm=mean(deg_brtzr,4);
deg_brtzrmz=zeros(size(deg_brtzrm));
for jj=1:size(deg_brtzrm,2)%bands
    for kk=1:size(deg_brtzrm,3)%time
         windeg=deg_brtzrm(:,jj,kk);
         deg_brtzrmz(:,jj,kk)=(deg_brtzrm(:,jj,kk)-mean(windeg(:)))/std(windeg(:));
    end
end
deg_brtzrmz(isnan(deg_brtzrmz))=0;
deg_brtzrm(isnan(deg_brtzrm))=0;
%deg_brtzrm(deg_final==0)=[];
percheckp(percnt)=max(deg_brtzrmz(:));
percheckn(percnt)=min(deg_brtzrmz(:));
end


%%%%real one
deg_brt=DEG_br;
deg_brtz=zeros(size(deg_brt));

deg_brtr=deg_brt;%(:,:,2:2:10,:);%zeros(size(deg_brt,1),size(deg_brt,2),10,size(deg_brt,4));

deg_brtzr=zeros(size(deg_brtr));

for jj=1:size(deg_brtr,2) %bands
    for ii=1:size(deg_brtr,3)%time
    
    for kk=1:size(deg_brtr,4)%subjs
        windeg=deg_brtr(:,jj,ii,kk);
        deg_brtzr(:,jj,ii,kk)=(deg_brtr(:,jj,ii,kk)-mean(windeg(:)))/std(windeg(:));
%         windeg1=deg_brr(:,jj,:,kk);
%         deg_brrz(:,jj,:,kk)=(deg_brr(:,jj,:,kk)-mean(windeg1(:)))/std(windeg1(:));
    end
    end
end

deg_brtzrm=mean(deg_brtzr,4);
deg_brtzrmz=zeros(size(deg_brtzrm));
for jj=1:size(deg_brtzrm,2)%bands
    for kk=1:size(deg_brtzrm,3)%time
         windeg=deg_brtzrm(:,jj,kk);
         deg_brtzrmz(:,jj,kk)=(deg_brtzrm(:,jj,kk)-mean(windeg(:)))/std(windeg(:));
    end
end
deg_brtzrmz(isnan(deg_brtzrmz))=0;
deg_sigp=zeros(size(deg_brtzrmz));
deg_sign=zeros(size(deg_brtzrmz));

for ii=1:size(deg_sigp,1)
    for jj=1:size(deg_sigp,2)
        for kk=1:size(deg_sigp,3)
    deg_sigp(ii,jj,kk)=1-sum(deg_brtzrmz(ii,jj,kk)>percheckp)/length(deg_brtzrmz(ii,jj,kk)>percheckp);
        deg_sign(ii,jj,kk)=1-sum(deg_brtzrmz(ii,jj,kk)<percheckn)/length(deg_brtzrmz(ii,jj,kk)>percheckn);

        end
    end
end
% deg_sigp(deg_final==0)=1;
% deg_sign(deg_final==0)=1;
