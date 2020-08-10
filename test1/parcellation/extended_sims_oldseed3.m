clear
clc
close all
%%aparc 
%SNR3-chunck1
cd('/imaging/rf02/parcellation/extended_simulations/singlesnr/sims70_offset0_seeds3_conns1_snr3_ext25/aparc')
[FileName,PathName,FilterIndex] = uigetfile('*.mat','Select the Mat-file','MultiSelect','on');
nsims=35;%length(FileName);
nnodes=68;
nsubs=17;
connmat=zeros(nnodes,nnodes,10,nsubs,nsims);%leakageSNR3,noiseSNR3,ideal,leakageSNR1,noiseSNR1 (COH 1st, imCOH2nd)
R0idn=zeros(nsubs,nsims);
R0idl=zeros(nsubs,nsims);
R0nl=zeros(nsubs,nsims);
% connmat=zeros(nnodes,nnodes,6,17,nsims);
for fcnt=1:14%nsims
a=load([PathName,FileName{fcnt}]);
connmat(:,:,[1,2,3,6,7,8],:,fcnt)=a.parc0_connmat(:,:,[1,2,3,6,7,8],:);%(:,:,:,:,fcnt)
end
%SNR1-chunck1
cd('/imaging/rf02/parcellation/extended_simulations/singlesnr/sims70_offset0_seeds3_conns1_snr1_ext25/aparc')
[FileName,PathName,FilterIndex] = uigetfile('*.mat','Select the Mat-file','MultiSelect','on');
for fcnt=1:14%nsims
a=load([PathName,FileName{fcnt}]);
connmat(:,:,[4,5,9,10],:,fcnt)=a.parc0_connmat(:,:,[1,2,6,7],:);%(:,:,:,:,fcnt)
end

%both SNRs
cd('/imaging/rf02/parcellation/extended_simulations/sims70_offset0_seeds3_conns1_bothsnrs_ext25/aparc')
[FileName,PathName,FilterIndex] = uigetfile('*.mat','Select the Mat-file','MultiSelect','on');
cc=0;
for fcnt=15:35%nsims
    cc=cc+1;
a=load([PathName,FileName{cc}]);
connmat(:,:,:,:,fcnt)=a.parc0_connmat;%(:,:,:,:,fcnt)
end

save('connmat_35sims.mat','connmat')

%%aparc_mod 
%SNR3-chunck1
cd('/imaging/rf02/parcellation/extended_simulations/singlesnr/sims70_offset0_seeds3_conns1_snr3_ext25/aparc_mod')
[FileName,PathName,FilterIndex] = uigetfile('*.mat','Select the Mat-file','MultiSelect','on');
nsims=35;%length(FileName);
nnodes=74;
nsubs=17;
connmat=zeros(nnodes,nnodes,10,nsubs,nsims);%leakageSNR3,noiseSNR3,ideal,leakageSNR1,noiseSNR1 (COH 1st, imCOH2nd)
R0idn=zeros(nsubs,nsims);
R0idl=zeros(nsubs,nsims);
R0nl=zeros(nsubs,nsims);
% connmat=zeros(nnodes,nnodes,6,17,nsims);
for fcnt=1:14%nsims
a=load([PathName,FileName{fcnt}]);
connmat(:,:,[1,2,3,6,7,8],:,fcnt)=a.parc1_connmat(:,:,[1,2,3,6,7,8],:);%(:,:,:,:,fcnt)
end
%SNR1-chunck1
cd('/imaging/rf02/parcellation/extended_simulations/singlesnr/sims70_offset0_seeds3_conns1_snr1_ext25/aparc_mod')
[FileName,PathName,FilterIndex] = uigetfile('*.mat','Select the Mat-file','MultiSelect','on');
for fcnt=1:14%nsims
a=load([PathName,FileName{fcnt}]);
connmat(:,:,[4,5,9,10],:,fcnt)=a.parc1_connmat(:,:,[1,2,6,7],:);%(:,:,:,:,fcnt)
end

%both SNRs
cd('/imaging/rf02/parcellation/extended_simulations/sims70_offset0_seeds3_conns1_bothsnrs_ext25/aparc_mod')
[FileName,PathName,FilterIndex] = uigetfile('*.mat','Select the Mat-file','MultiSelect','on');
cc=0;
for fcnt=15:35%nsims
    cc=cc+1;
a=load([PathName,FileName{cc}]);
connmat(:,:,:,:,fcnt)=a.parc1_connmat;%(:,:,:,:,fcnt)
end

save('connmat_35sims.mat','connmat')

%%aparc2009 
%SNR3-chunck1
cd('/imaging/rf02/parcellation/extended_simulations/singlesnr/sims70_offset0_seeds3_conns1_snr3_ext25/aparc2009')
[FileName,PathName,FilterIndex] = uigetfile('*.mat','Select the Mat-file','MultiSelect','on');
nsims=35;%length(FileName);
nnodes=148;
nsubs=17;
connmat=zeros(nnodes,nnodes,10,nsubs,nsims);%leakageSNR3,noiseSNR3,ideal,leakageSNR1,noiseSNR1 (COH 1st, imCOH2nd)
R0idn=zeros(nsubs,nsims);
R0idl=zeros(nsubs,nsims);
R0nl=zeros(nsubs,nsims);
% connmat=zeros(nnodes,nnodes,6,17,nsims);
for fcnt=1:14%nsims
a=load([PathName,FileName{fcnt}]);
connmat(:,:,[1,2,3,6,7,8],:,fcnt)=a.parc2_connmat(:,:,[1,2,3,6,7,8],:);%(:,:,:,:,fcnt)
end
%SNR1-chunck1
cd('/imaging/rf02/parcellation/extended_simulations/singlesnr/sims70_offset0_seeds3_conns1_snr1_ext25/aparc2009')
[FileName,PathName,FilterIndex] = uigetfile('*.mat','Select the Mat-file','MultiSelect','on');
for fcnt=1:14%nsims
a=load([PathName,FileName{fcnt}]);
connmat(:,:,[4,5,9,10],:,fcnt)=a.parc2_connmat(:,:,[1,2,6,7],:);%(:,:,:,:,fcnt)
end

%both SNRs
cd('/imaging/rf02/parcellation/extended_simulations/sims70_offset0_seeds3_conns1_bothsnrs_ext25/aparc2009')
[FileName,PathName,FilterIndex] = uigetfile('*.mat','Select the Mat-file','MultiSelect','on');
cc=0;
for fcnt=15:35%nsims
    cc=cc+1;
a=load([PathName,FileName{cc}]);
connmat(:,:,:,:,fcnt)=a.parc2_connmat;%(:,:,:,:,fcnt)
end

save('connmat_35sims.mat','connmat')

%%aparc2009_mod
%SNR3-chunck1
cd('/imaging/rf02/parcellation/extended_simulations/singlesnr/sims70_offset0_seeds3_conns1_snr3_ext25/aparc2009_mod')
[FileName,PathName,FilterIndex] = uigetfile('*.mat','Select the Mat-file','MultiSelect','on');
nsims=35;%length(FileName);
nnodes=74;
nsubs=17;
connmat=zeros(nnodes,nnodes,10,nsubs,nsims);%leakageSNR3,noiseSNR3,ideal,leakageSNR1,noiseSNR1 (COH 1st, imCOH2nd)
R0idn=zeros(nsubs,nsims);
R0idl=zeros(nsubs,nsims);
R0nl=zeros(nsubs,nsims);
% connmat=zeros(nnodes,nnodes,6,17,nsims);
for fcnt=1:14%nsims
a=load([PathName,FileName{fcnt}]);
connmat(:,:,[1,2,3,6,7,8],:,fcnt)=a.parc3_connmat(:,:,[1,2,3,6,7,8],:);%(:,:,:,:,fcnt)
end
%SNR1-chunck1
cd('/imaging/rf02/parcellation/extended_simulations/singlesnr/sims70_offset0_seeds3_conns1_snr1_ext25/aparc2009_mod')
[FileName,PathName,FilterIndex] = uigetfile('*.mat','Select the Mat-file','MultiSelect','on');
for fcnt=1:14%nsims
a=load([PathName,FileName{fcnt}]);
connmat(:,:,[4,5,9,10],:,fcnt)=a.parc3_connmat(:,:,[1,2,6,7],:);%(:,:,:,:,fcnt)
end

%both SNRs
cd('/imaging/rf02/parcellation/extended_simulations/sims70_offset0_seeds3_conns1_bothsnrs_ext25/aparc2009_mod')
[FileName,PathName,FilterIndex] = uigetfile('*.mat','Select the Mat-file','MultiSelect','on');
cc=0;
for fcnt=15:35%nsims
    cc=cc+1;
a=load([PathName,FileName{cc}]);
connmat(:,:,:,:,fcnt)=a.parc3_connmat;%(:,:,:,:,fcnt)
end

save('connmat_35sims.mat','connmat')

%%RG
%SNR3-chunck1
cd('/imaging/rf02/parcellation/extended_simulations/singlesnr/sims70_offset0_seeds3_conns1_snr3_ext25/RG')
[FileName,PathName,FilterIndex] = uigetfile('*.mat','Select the Mat-file','MultiSelect','on');
nsims=35;%length(FileName);
nnodes=70;
nsubs=17;
connmat=zeros(nnodes,nnodes,10,nsubs,nsims);%leakageSNR3,noiseSNR3,ideal,leakageSNR1,noiseSNR1 (COH 1st, imCOH2nd)
R0idn=zeros(nsubs,nsims);
R0idl=zeros(nsubs,nsims);
R0nl=zeros(nsubs,nsims);
% connmat=zeros(nnodes,nnodes,6,17,nsims);
for fcnt=1:14%nsims
a=load([PathName,FileName{fcnt}]);
connmat(:,:,[1,2,3,6,7,8],:,fcnt)=a.parc4_connmat(:,:,[1,2,3,6,7,8],:);%(:,:,:,:,fcnt)
end
%SNR1-chunck1
cd('/imaging/rf02/parcellation/extended_simulations/singlesnr/sims70_offset0_seeds3_conns1_snr1_ext25/RG')
[FileName,PathName,FilterIndex] = uigetfile('*.mat','Select the Mat-file','MultiSelect','on');
for fcnt=1:14%nsims
a=load([PathName,FileName{fcnt}]);
connmat(:,:,[4,5,9,10],:,fcnt)=a.parc4_connmat(:,:,[1,2,6,7],:);%(:,:,:,:,fcnt)
end

%both SNRs
cd('/imaging/rf02/parcellation/extended_simulations/sims70_offset0_seeds3_conns1_bothsnrs_ext25/RG')
[FileName,PathName,FilterIndex] = uigetfile('*.mat','Select the Mat-file','MultiSelect','on');
cc=0;
for fcnt=15:35%nsims
    cc=cc+1;
a=load([PathName,FileName{cc}]);
connmat(:,:,:,:,fcnt)=a.parc4_connmat;%(:,:,:,:,fcnt)
end

save('connmat_35sims.mat','connmat')

% for subcnt=1:nsubs
% pm0=connmat(:,:,:,subcnt);%squeeze(mean(connmat,4));
% pm0id=((pm0(:,:,3)+pm0(:,:,3)')/2)+2*eye(size(pm0,1),size(pm0,2));
% pm0n=((pm0(:,:,2)+pm0(:,:,2)')/2)+2*eye(size(pm0,1),size(pm0,2));
% pm0l=((pm0(:,:,1)+pm0(:,:,1)')/2)+2*eye(size(pm0,1),size(pm0,2));
% pm0id(pm0id==2)=[];
% pm0id=reshape(pm0id,nnodes-1,nnodes);
% pm0n(pm0n==2)=[];
% pm0n=reshape(pm0n,nnodes-1,nnodes);
% pm0l(pm0l==2)=[];
% pm0l=reshape(pm0l,nnodes-1,nnodes);
% R0idn(subcnt,fcnt)=corr2(pm0id,pm0n);
% R0idl(subcnt,fcnt)=corr2(pm0id,pm0l);
% R0nl(subcnt,fcnt)=corr2(pm0n,pm0l);
% end


%%aparc mod
% % % [FileName,PathName,FilterIndex] = uigetfile('*.mat','Select the Mat-file','MultiSelect','on');
% % % nsims=length(FileName);
% % % nsubs=17;
% % % nnodes=74;
% % % R1idn=zeros(nsubs,nsims);
% % % R1idl=zeros(nsubs,nsims);
% % % R1nl=zeros(nsubs,nsims);
% % % % connmat=zeros(nnodes,nnodes,6,17,nsims);
% % % for fcnt=1:nsims
% % % a=load([PathName,FileName{fcnt}]);
% % % connmat=a.parc1_connmat;%(:,:,:,:,fcnt)
% % % for subcnt=1:nsubs
% % % pm1=connmat(:,:,:,subcnt);%squeeze(mean(connmat,4));
% % % pm1id=((pm1(:,:,3)+pm1(:,:,3)')/2)+2*eye(size(pm1,1),size(pm1,2));
% % % pm1n=((pm1(:,:,2)+pm1(:,:,2)')/2)+2*eye(size(pm1,1),size(pm1,2));
% % % pm1l=((pm1(:,:,1)+pm1(:,:,1)')/2)+2*eye(size(pm1,1),size(pm1,2));
% % % pm1id(pm1id==2)=[];
% % % pm1id=reshape(pm1id,nnodes-1,nnodes);
% % % pm1n(pm1n==2)=[];
% % % pm1n=reshape(pm1n,nnodes-1,nnodes);
% % % pm1l(pm1l==2)=[];
% % % pm1l=reshape(pm1l,nnodes-1,nnodes);
% % % R1idn(subcnt,fcnt)=corr2(pm1id,pm1n);
% % % R1idl(subcnt,fcnt)=corr2(pm1id,pm1l);
% % % R1nl(subcnt,fcnt)=corr2(pm1n,pm1l);
% % % end
% % % end