clear
clc
close all
% addpath('/home/rf02/rezvan/test1/BCT/2015_01_25 BCT')
%in_path='C:\Users\rf02\Documents\Rezvan\PhDproject\semnet\rezvan\rezvan\rezvan_meg_final_meglab007\data\';%'/imaging/rf02/TypLexMEG/';
in_path='/home/rf02/PC_backup/PhDproject/semnet/rezvan/rezvan/rezvan_meg_final_meglab007/data/';
load('/imaging/rf02/semnet_git/SemNet/semnet4semloc/word_concrete_matchedca50_2.mat')
load('/imaging/rf02/semnet_git/SemNet/semnet4semloc/word_eabs_matchedca50_2.mat')

list_all = {'MEG_0032_MEG16002_29_2_18_11_1_LDblock.mat',
    'MEG0034_MEG16002_41_2_19_13_57_LDblock.mat',
    'MEG16_0035_2_MEG16002_26_2_22_17_46_LDblock.mat',
    'MEG16_0042_MEG16002_24_2_29_11_53_LDblock.mat',
    'MEG16_0045_MEG16002_27_3_3_13_46_LDblock.mat',
    'MEG16_0052_2_MEG16002_24_3_10_13_59_LDblock.mat',
    'MEG16_0056_MEG16002_36_3_14_15_36_LDblock.mat',
    'MEG16_0069_MEG16002_30_4_5_15_51_LDblock.mat',
    'MEG16_0070_MEG16002_28_4_7_9_17_LDblock.mat',
    'MEG16_0072_MEG16002_23_4_8_15_22_LDblock.mat',
    'MEG16_0073_MEG16002_26_4_11_9_24_LDblock.mat',
    'MEG16_0075_MEG16002_33_4_11_15_57_LDblock.mat',
    'MEG16_0078_MEG16002_23_4_14_11_23_LDblock.mat',
    'MEG16_0082_MEG16002_22_4_18_11_32_LDblock.mat',
    'MEG16_0086_MEG16002_29_4_22_15_43_LDblock.mat',
    'MEG16_0097_MEG16002_23_5_12_11_41_LDblock.mat',
    'MEG16_0122_MEG16002_22_7_7_9_27_LDblock.mat',
    %'Meg16_0123_MEG16002_23_7_8_9_41_LDblock.mat',
    'MEG16_0125_MEG16002_24_7_12_15_48_LDblock.mat',
    };


n_subj=length(list_all);
rtmat=zeros(n_subj,8);
crrate=zeros(n_subj,8);
checkbads=zeros(n_subj,3,500);
p=[];
for cnt=1:n_subj
    cnt
    fname=[in_path,list_all{cnt}];
    load(fname)
    [aa,bb]=sort(wordlist);
    respmats=respmat(bb,:);
    checkbads(cnt,:,:)=respmats(:,[3,6,7])';

end
cb=squeeze(mean(checkbads,1))';
goods=cb(:,1)>0.5;
bads=cb(:,1)<=0.5;
auds=cb(:,2)==2;
audsl=find(cb(:,2)==2);

audcr=cb(auds,1);
[audaa,audbb]=sort(audcr);
audsls=audsl(audbb);
audws=aa(auds);
audwsaa=audws(audbb);
ideals=cb(:,1)==1;

c=cb(bads,:);
visbads=sum(c(:,2)==1);
audbads=sum(c(:,2)==2);
hndbads=sum(c(:,2)==3);
pwbads=sum(c(:,2)==6);
% idealsvis=ideals;
% idealsvis(cb(:,2)~=1)=0;
% wiv=find(idealsvis==1);
% wiv=wiv(1:audbads-visbads);
% goods(wiv)=0;
% 
% idealshnd=ideals;
% idealshnd(cb(:,2)~=3)=0;
% wih=find(idealshnd==1);
% wih=wih(1:audbads-hndbads);
% goods(wih)=0;
% 
% idealspw=ideals;
% idealspw(cb(:,2)~=6)=0;
% wip=find(idealspw==1);
% wip=wip(1:audbads-pwbads);
% goods(wip)=0;
% 
% ci=cb(ideals,:);
% aa(ans);
audini=audbads;


% while min(per(:))<0.05
%     audbads=audbads+1;
%     goods(audsls(audbads))=0;
for cnt=1:n_subj
    cnt
    
    fname=[in_path,list_all{cnt}];
    load(fname)
    con_ind=zeros(length(wordlist),1);
    abs_ind=zeros(length(wordlist),1);
    for con_cnt=1:length(con_list)
        thiscon=find(strcmp(wordlist,con_list(con_cnt)));
        con_ind(thiscon)=1;
    end
    for abs_cnt=1:length(eabs_list)
        thisabs=find(strcmp(wordlist,eabs_list(abs_cnt)));
        abs_ind(thisabs)=1;
    end
    
    respmat(con_ind==1,6)=20;%concrete
    respmat(abs_ind==1,6)=21;%emot abstract
    
    [aa,bb]=sort(wordlist);
    respmats=respmat(bb,:);
    respmatsg=respmats(goods==1,:);
%     checkbads(cnt,:,:)=respmatsg(:,[3,6])';
    respmat_cr=respmats(respmats(:,3)==1,:);
    
%     vismat=respmat_cr(respmat_cr(:,6)==1,:);rtmat(cnt,1)=mean(vismat(:,2));crrate(cnt,1)=size(vismat,1)/(50-visbads);%50;%(50-audbads);%
%     audmat=respmat_cr(respmat_cr(:,6)==2,:);rtmat(cnt,2)=mean(audmat(:,2));crrate(cnt,2)=size(audmat,1)/(50-audbads);%(50-audbads);%(50-sum(c(:,2)==2));
%     hndmat=respmat_cr(respmat_cr(:,6)==3,:);rtmat(cnt,3)=mean(hndmat(:,2));crrate(cnt,3)=size(hndmat,1)/(50-hndbads);%(50-audbads);%(50-sum(c(:,2)==3));
%     emotmat=respmat_cr(respmat_cr(:,6)==4,:);rtmat(cnt,4)=mean(emotmat(:,2));crrate(cnt,4)=size(emotmat,1)/(50-sum(c(:,2)==4));%
%     neutmat=respmat_cr(respmat_cr(:,6)==5,:);rtmat(cnt,5)=mean(neutmat(:,2));crrate(cnt,5)=size(neutmat,1)/(50-sum(c(:,2)==5));%(50-sum(c(:,2)==5));
%     pwmat=respmat_cr(respmat_cr(:,6)==6,:);rtmat(cnt,6)=mean(pwmat(:,2));crrate(cnt,6)=size(pwmat,1)/(50-pwbads);%(50-audbads);%(50-sum(c(:,2)==6));
%     
    conmat=respmat_cr(respmat_cr(:,6)==20,:);rtmat(cnt,7)=mean(conmat(:,2));crrate(cnt,7)=size(conmat,1)/50;%50;%(50-audbads);%
    absmat=respmat_cr(respmat_cr(:,6)==21,:);rtmat(cnt,8)=mean(absmat(:,2));crrate(cnt,8)=size(absmat,1)/50;%50;%(50-audbads);%

   
end
% [h,p(1)]=ttest(rtmat(:,1)-rtmat(:,2));%vis aud
% [h,p(2)]=ttest(rtmat(:,1)-rtmat(:,3));%vis hnd
% [h,p(3)]=ttest(rtmat(:,2)-rtmat(:,3));%aud hnd
% [h,p(4)]=ttest(rtmat(:,2)-rtmat(:,4));%aud hnd
rtmat=rtmat(:,[7,8]);
prt=eye(size(rtmat,2),size(rtmat,2));
for ii=1:size(rtmat,2)
    for jj=setdiff(1:size(rtmat,2),ii)
%         [h,prt(ii,jj)]=ttest(rtmat(:,ii)-rtmat(:,jj));
prt(ii,jj)=signrank(rtmat(:,ii)-rtmat(:,jj));
    end
end
crrate=crrate(:,[7,8]);
per=eye(size(crrate,2),size(crrate,2));

for ii=1:size(crrate,2)
    for jj=setdiff(1:size(crrate,2),ii)
%         [h,per(ii,jj)]=ttest(crrate(:,ii)-crrate(:,jj));
per(ii,jj)=signrank(crrate(:,ii)-crrate(:,jj));
    end
end
% end
% rtmat=rtmat(:,[1,2,3,6]);
boxplot(rtmat,'notch','on','labels',{'Concrete','Abstract'})%,'boxstyle','filled')

% boxplot(rtmat,'notch','on','labels',{'Visual','Auditory','Hand Action','Emotional Abs','Neutral Abs','Pseudowords'})%,'boxstyle','filled')

figure,
% boxplot(crrate,'notch','on','labels',{'Visual','Auditory','Hand Action','Emotional Abs','Neutral Abs','Pseudowords'})%,'boxstyle','filled')
boxplot(crrate,'notch','on','labels',{'Concrete','Abstract'})%,'boxstyle','filled')

