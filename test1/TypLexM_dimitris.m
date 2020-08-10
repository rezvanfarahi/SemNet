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

fname_a='label_mne_ctf_mn_ts_abstract.mat';
fname_c='label_mne_ctf_mn_ts_concrete.mat';

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
close all
for cnt=1:17
     cnt
    fname=[in_path,list_all{cnt},fname_a];
    load(fname)
    
    fname=[in_path,list_all{cnt},fname_c];
    load(fname)
%     

b=0;
for ii=1
tc=double(squeeze(label_mne_ctf_mn_ts_concrete(:,ii+b,750:950)));
ta=double(squeeze(label_mne_ctf_mn_ts_abstract(:,ii+b,750:950)));
ts=[tc;ta]';
C=corr(ts);
figure,
imagesc((C))
%colormap(hot)
colorbar
end
end