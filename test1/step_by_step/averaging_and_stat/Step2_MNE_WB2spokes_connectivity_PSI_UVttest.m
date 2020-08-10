clear
clc
close all
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

fname_a='mean_coh_abstract_multitaper_8regions.mat';
fname_c='mean_coh_concrete_multitaper_8regions.mat';
Matx_c=zeros(8,8,2,2,17);
Matx_a=zeros(8,8,2,2,17);

fname=[in_path,list_all{1},fname_a];
load(fname)
% fname_s='plv_subtract_cwtmorlet.mat';
n_subj=length(list_all);
node_perc=[0.05:0.05:0.75];

tr=zeros(n_subj,1);
trb=zeros(n_subj,1);
%labellist = ['lh.posterior_fusiform','rh.posterior_fusiform', 'lh.lateraloccipital','rh.lateraloccipital', 'lh.ATLsmall','rh.ATLsmall', 'lh.entorhinal', 'rh.entorhinal']#labellist = [ 'leftATL-lh','rightATL-rh']#, 'medialtempright-rh','medialtempleft-lh']'lh.postcentral','rh.postcentral',

load('/imaging/rf02/TypLexMEG/icaanalysis_results/jpg/time_resolved/ROI_connectivity/graph/labels_46.mat')
load('/imaging/rf02/TypLexMEG/icaanalysis_results/jpg/time_resolved/ROI_connectivity/graph/node_angles46.mat')
close all
for cnt=1:17
    cnt
    fname=[in_path,list_all{cnt},fname_c];
    load(fname)
    Matx_c(:,:,:,1,cnt)=mean(con_cncrt(:,:,:,650:850),4);
    Matx_c(:,:,:,2,cnt)=mean(con_cncrt(:,:,:,750:950),4);
    
    fname=[in_path,list_all{cnt},fname_a];
    load(fname)
    Matx_a(:,:,:,1,cnt)=mean(con_abs(:,:,:,650:850),4);
    Matx_a(:,:,:,2,cnt)=mean(con_abs(:,:,:,750:950),4);
    %
    
end
for ii=1:size(Matx_a,5)%subj
    for jj=1:size(Matx_a,4)%time
        for kk=1:size(Matx_a,3)%bands
            for mm=1:size(Matx_a,2)
                for nn=1:size(Matx_a,1)
                    if mm>nn
                        Matx_a(nn,mm,kk,jj,ii)=Matx_a(mm,nn,kk,jj,ii);
                        Matx_c(nn,mm,kk,jj,ii)=Matx_c(mm,nn,kk,jj,ii);
                    end
                end
            end
        end
    end
end
ttest_a=zeros(8,8,2,2);
ttest_c=zeros(8,8,2,2);
ttest_s=zeros(8,8,4,2);
for jj=1:size(ttest_a,4)%time
    for kk=1:size(ttest_a,3)%bands
        for mm=1:size(ttest_a,2)
            for nn=1:size(ttest_a,1)
                [h,o,ci,stats]=ttest(squeeze(Matx_a(nn,mm,kk,jj,:)));
                ttest_a(nn,mm,kk,jj)=stats.tstat;
                [h,o,ci,stats]=ttest(squeeze(Matx_c(nn,mm,kk,jj,:)));
                ttest_c(nn,mm,kk,jj)=stats.tstat;
                [h,o,ci,stats]=ttest(squeeze(Matx_c(nn,mm,kk,jj,:)-Matx_a(nn,mm,kk,jj,:)));
                ttest_s(nn,mm,kk,jj)=stats.tstat;
            end
        end
    end
end

ttest_a(isnan(ttest_a))=0;
ttest_c(isnan(ttest_c))=0;
Ma=Matx_a;
Mc=Matx_c;


Matx_a(Matx_a>0)=1;
Matx_a(Matx_a<0)=0;
ma=mean(Matx_a,5);
Matx_c(Matx_c>0)=1;
Matx_c(Matx_c<0)=0;
mc=mean(Matx_c,5);
for ii=1:size(ma,3)%bands
    for jj=1:size(ma,4)%time
        figure, imagesc(squeeze(ttest_s(5:end,1:4,ii,jj)),[-2,2]), colorbar
        %figure,imagesc(squeeze(ttest_a(5:end,1:4,ii,jj)),[-2,2]), colorbar
    end
end

