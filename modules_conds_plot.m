clear
clc
close all
addpath('/imaging/rf02/Semnet/MIToolbox_master/matlab')
addpath('/imaging/rf02/Semnet/GenLouvain_master')
addpath('/imaging/rf02/Semnet/BCT/2017_01_15_BCT')
out_path='/imaging/rf02/Semnet/graph_analysis/';
data_path='/imaging/rf02/Semnet';
list_all =  {'/meg16_0030/160216/', 
            '/meg16_0032/160218/', 
            '/meg16_0034/160219/', 
            '/meg16_0035/160222/',
            '/meg16_0042/160229/', 
            '/meg16_0045/160303/', 
            '/meg16_0052/160310/', 
            '/meg16_0056/160314/',
            '/meg16_0069/160405/', 
            '/meg16_0070/160407/', 
            '/meg16_0072/160408/', 
            '/meg16_0073/160411/', 
            '/meg16_0075/160411/', 
            '/meg16_0078/160414/', 
            '/meg16_0082/160418/', 
            '/meg16_0086/160422/', 
            '/meg16_0097/160512/',  
            '/meg16_0122/160707/', 
            '/meg16_0125/160712/', 
            };

% fname_vis='pacel_connmat_Visual_aparc_mod.mat';
% fname_hnd='pacel_connmat_Hand_aparc_mod.mat';
% fname_aud='pacel_connmat_Hear_aparc_mod.mat';
fname_conds={'pacel_connmat_maxctf_Visual_aparc_mod.mat','pacel_connmat_maxctf_Hear_aparc_mod.mat','pacel_connmat_maxctf_Hand_aparc_mod.mat'};
nsubs=length(list_all);
nrands=100;
tau=0.1;%nrands/10;
reps=100;
nnodes=74;
timesmin={'50ms','150ms','250ms'};
timesmax={'250ms','350ms','450ms'};
ntimes=length(timesmin);
nbands=4;
nconds=length(fname_conds);


% gamma=2.6;
% hnd_connmat=zeros(nnodes,nnodes,2,nsubs);
grand_connmat=zeros(nsubs,nnodes,nnodes,nbands,ntimes,nconds);
grand_connmat_thresh=zeros(nsubs,nnodes,nnodes,nbands,ntimes,nconds);
grand_parmat=zeros(nsubs,nnodes,nbands,ntimes,nconds);
out_file_modules=[out_path,'wmodules_mi_maxctf_3conds_up4_3bands.mat'];
load(out_file_modules)
out_file_modules=[out_path,'finalmodules_mi_maxctf_3conds_up4_3bands.mat'];
load(out_file_modules)

modunique=unique(finalmodules);
nmodules=length(modunique);
Cof=0.25;

for ii=1:nsubs
    ii
connmat1=zeros(nnodes,nnodes,nbands,ntimes);
connmat_thresh=zeros(nnodes,nnodes,nbands,ntimes);
parmat=zeros(nnodes,nbands,ntimes);

for condcnt=1:nconds
    condcnt
    fnamecond_full=[data_path,list_all{ii},fname_conds{condcnt}];
    load(fnamecond_full)
    
    for jj=1:size(connmat,3)%freqs
        for kk=1:size(connmat,4)%times
            aa=squeeze(connmat(:,:,jj,kk));
            aa1=aa+aa'+eye(size(aa));
            aa2=aa+aa';
%             aa(zscore(aa)
% bb=tril(aa,-1);
% cc=bb;
% cc(cc>0)=zscore(bb(bb>0));
% cc(cc<1.96)=0;
            connmat1(:,:,jj,kk)=aa2(wmodules,wmodules);
%             connmat_thresh(:,:,jj,kk)=threshold_proportional(aa2(wmodules,wmodules),0.15);
%             parmat(:,jj,kk)=participation_coef(aa2,finalmodules);
        end
    end
    plrb=connmat1(:);
    [aab,bbb]=sort(plrb,'descend'); 
    aab2=aab(aab>0);
    trb=aab2(round(Cof*length(aab2)));
    for jj=1:size(connmat,3)%freqs
        for kk=1:size(connmat,4)%times
            aa=squeeze(connmat(:,:,jj,kk));
            aa1=aa+aa'+eye(size(aa));
            aa2=aa+aa';
%             aa(zscore(aa)
bb=tril(aa2,-1);
cc=bb;
cc(cc>0)=zscore(bb(bb>0));
cc=cc+cc';
% cc(cc<1.96)=0;
            connmat1(:,:,jj,kk)=aa1(wmodules,wmodules);
            connmat_thresh(:,:,jj,kk)=aa2(wmodules,wmodules);%aa2(wmodules,wmodules);%threshold_absolute(aa2(wmodules,wmodules),trb);
            parmat(:,jj,kk)=module_degree_zscore(squeeze(connmat_thresh(:,:,jj,kk)),finalmodules);%participation_coef(squeeze(connmat_thresh(:,:,jj,kk)),finalmodules);%gateway_coef_sign(cc(wmodules,wmodules),finalmodules,1);%participation_coef(aa2(wmodules,wmodules),finalmodules);%
        end
    end
    
    grand_connmat(ii,:,:,:,:,condcnt)=connmat1;
%     connmat_thresh(connmat_thresh>0)=1;
    grand_connmat_thresh(ii,:,:,:,:,condcnt)=connmat_thresh;
    grand_parmat(ii,:,:,:,condcnt)=parmat;
    
end
end


%% create module connmat
thresh_connmat_modules=zeros(nsubs,nmodules,nmodules,nbands,ntimes,nconds);
betweenness_modules=zeros(nsubs,nmodules,nbands,ntimes,nconds);

for ii=1:nsubs
    for jj=1:nbands
        for kk=1:ntimes
            for ll=1:nconds
                thisM=squeeze(grand_connmat_thresh(ii,:,:,jj,kk,ll));         
                for mm=1:nmodules
                    for nn=1:nmodules
                        thismod=thisM(finalmodules==mm,finalmodules==nn);
                        if mm==nn
                        thresh_connmat_modules(ii,mm,nn,jj,kk,ll)=sum(thismod(:))/(size(thismod,1)*(size(thismod,1)-1));
                        else
                            thresh_connmat_modules(ii,mm,nn,jj,kk,ll)=mean(thismod(:));
                        end
                    end
                end
                betweenness_modules(ii,:,jj,kk,ll)=strengths_und(squeeze(thresh_connmat_modules(ii,:,:,jj,kk,ll)));%betweenness_wei(squeeze(thresh_connmat_modules(ii,:,:,jj,kk,ll)));
            end
        end
    end
end

                        

% out_file_grandmat=[out_path,'grand_mat_3conds.mat'];
% load(out_file_grandmat)
% out_file_Dgrand=[out_path,'D_grand_3conds.mat'];
% load(out_file_Dgrand)
% out_file_ciugrand=[out_path,'ciu_grand_3conds.mat'];
% load(out_file_ciugrand)

mean_connmat=squeeze(mean(grand_connmat_thresh,1));
std_connmat=squeeze(std(thresh_connmat_modules,1));

mconnmat_alpha=squeeze(mean_connmat(:,:,2,:,:));
mconnmat_gamma=squeeze(mean_connmat(:,:,4,:,:));
sconnmat_alpha=squeeze(std_connmat(:,:,2,:,:));
sconnmat_gamma=squeeze(std_connmat(:,:,4,:,:));

thiscond=[2,3,3];
refcond=[1,1,2];
figure,
kk=0;
for ii=1:nconds
    for jj=1:ntimes
        kk=kk+1;
        subplot(ntimes,nconds,kk)
        imagesc(squeeze(mconnmat_alpha(:,:,jj,refcond(ii)))-squeeze(mconnmat_alpha(:,:,jj,thiscond(ii))),[-0.01 0.01])
    end
end
% colormap(jet)
% colorbar

figure,
kk=0;
for ii=1:nconds
    for jj=1:ntimes
        kk=kk+1;
        subplot(ntimes,nconds,kk)
        bar(squeeze(sum(mconnmat_alpha(:,:,jj,refcond(ii)),2))-squeeze(sum(mconnmat_alpha(:,:,jj,thiscond(ii)),2)))%,[-0.01 0.01])
    end
end




% figure,=1:nconds
% kk=0;
% thiscond1=[2,3;1,2;1,3];
% refcond1=[1,2,3];
% for ii=1:nconds
%     for jj=1:ntimes
%         kk=kk+1;
%         subplot(ntimes,nconds,kk)
%         imagesc(squeeze(mconnmat_alpha(:,:,jj,refcond1(ii)))-squeeze(mean(mconnmat_alpha(:,:,jj,thiscond1(ii,:)),4)),[-0.01 0.01])
%     end
% end
% colormap(jet)
% colorbar





figure,
kk=0;
for ii=1:nconds
    for jj=1:ntimes
        kk=kk+1;
        subplot(ntimes,nconds,kk)
        imagesc(squeeze(mconnmat_gamma(:,:,jj,refcond(ii)))-squeeze(mconnmat_gamma(:,:,jj,thiscond(ii))),[-0.01 0.01])
    end
end
% colormap(jet)
% colorbar

% figure,
% kk=0;
% thiscond1=[2,3;1,2;1,3];
% refcond1=[1,2,3];
% for ii=1:nconds
%     for jj=1:ntimes
%         kk=kk+1;
%         subplot(ntimes,nconds,kk)
%         imagesc(squeeze(mconnmat_gamma(:,:,jj,refcond1(ii)))-squeeze(mean(mconnmat_gamma(:,:,jj,thiscond1(ii,:)),4)),[-0.01 0.01])
%     end
% end

figure,
kk=0;
for ii=1:nconds
    for jj=1:ntimes
        kk=kk+1;
        subplot(ntimes,nconds,kk)
        bar(squeeze(sum(mconnmat_gamma(:,:,jj,refcond(ii)),2))-squeeze(sum(mconnmat_gamma(:,:,jj,thiscond(ii)),2)))%,[-0.01 0.01])
        
    end
end

figure,
kk=0;
for ii=1:nconds
    for jj=1:ntimes
        kk=kk+1;
        subplot(ntimes,nconds,kk)
        imagesc(squeeze(mconnmat_alpha(:,:,jj,ii)))
    end
end
colormap(hot)
        
figure,
kk=0;
for ii=1:nconds
    for jj=1:ntimes
        kk=kk+1;
        subplot(ntimes,nconds,kk)
        imagesc(squeeze(mconnmat_gamma(:,:,jj,ii)))
    end
end       
colormap(hot)
        
load('/imaging/rf02/Semnet/graph_analysis/labels_aparc_mod.mat')
lmod=label_names(wmodules,:);
% out_file_grandconnmat=[out_path,'grand_connmat_mi_maxctf_modules3conds_3bands.mat'];
% save(out_file_grandconnmat,'grand_connmat')
% out_file_grandconnmat=[out_path,'grand_connmat_mi_maxctf_thresh_modules3conds_3bands.mat'];
% save(out_file_grandconnmat,'grand_connmat_thresh')
% out_file_grandconnmat=[out_path,'grand_participation_mi_maxctf_nodes_modules3conds_3bands.mat'];
% save(out_file_grandconnmat,'grand_parmat')
% out_file_grandconnmat=[out_path,'grand_modconnmat_mi_maxctf_thresh_modules3conds_3bands.mat'];
% save(out_file_grandconnmat,'thresh_connmat_modules')

 
% thisfreq=3;
% thisev=3;
% refev=1;
% pval_mat=zeros(size(thresh_connmat_modules,2),size(thresh_connmat_modules,3),size(thresh_connmat_modules,5));
% for ii=1:size(pval_mat,1)
%     for jj=1:size(pval_mat,2)
%         for kk=1:size(pval_mat,3)
%             pval_mat(ii,jj,kk)=signrank(squeeze(thresh_connmat_modules(:,ii,jj,thisfreq,kk,refev))-squeeze(thresh_connmat_modules(:,ii,jj,thisfreq,kk,thisev)));
%         end
%     end
% end
% pval_mat(pval_mat>0.05)=1;

% thisfreq=[];
thisev=[2,3,3];
refev=[1,1,2];
pval_mat=zeros(nmodules,ntimes,4,nconds);
for ii=1:size(pval_mat,1)
        for jj=1:size(pval_mat,2)
            for kk=1:size(pval_mat,3)
                for ll=1:size(pval_mat,4)
            pval_mat(ii,jj,kk,ll)=signrank(squeeze(betweenness_modules(:,ii,kk,jj,refev(ll)))-squeeze(betweenness_modules(:,ii,kk,jj,thisev(ll))));
                end
            end
        end
end
pval_mat(pval_mat>0.05)=1;

FDR=[];
thisev=2;
refev=3;
for thisfreq=1:4;
pval_mat_par=zeros(nnodes,ntimes);
for ii=1:size(pval_mat_par,1)
    for jj=1:size(pval_mat_par,2)
            pval_mat_par(ii,jj)=signrank(squeeze(grand_parmat(:,ii,thisfreq,jj,refev))-squeeze(grand_parmat(:,ii,thisfreq,jj,thisev)));
    end
end
for ii=1:ntimes
FDR(:,ii,thisfreq) = mafdr(pval_mat_par(:,ii));
end
FDR(FDR>=0.05)=1;
min(FDR(:))
end
% pval_mat_par(pval_mat_par>0.05)=1;
            
            
            
            


