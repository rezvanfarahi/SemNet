clear
clc
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
            '/meg16_0123/160708/', 
            '/meg16_0125/160712/', 
            };

fname_vis='pacel_connmat_Visual_aparc_mod.mat';
fname_hnd='pacel_connmat_Hand_aparc_mod.mat';

nsubs=length(list_all);
nrands=100;
tau=0.1;%nrands/10;
reps=100;
nnodes=74;
timesmin={'50ms','150ms','250ms'};
timesmax={'250ms','350ms','450ms'};
ntimes=length(timesmin);
nbands=2;
nconds=2;


% gamma=2.6;
% hnd_connmat=zeros(nnodes,nnodes,2,nsubs);
ngammas=length(1:0.1:3);
grand_mat_sub_time_band_cond=zeros(nnodes,nsubs*ntimes*nbands*nconds,ngammas);
D_grand=zeros(nnodes,nnodes,ngammas);
ciu_grand=zeros(nnodes,ngammas);
gcnt=0;
for gamma=1:0.1:3
    gcnt=gcnt+1;
    gamma
    ciu_vis_alpha=zeros(nnodes,nsubs,ntimes);
ciu_vis_gamma=zeros(nnodes,nsubs,ntimes);
ciu_hnd_alpha=zeros(nnodes,nsubs,ntimes);
ciu_hnd_gamma=zeros(nnodes,nsubs,ntimes);
for thist=1:ntimes
    thist
    vis_alpha_connmat=zeros(nnodes,nnodes,nsubs);
vis_gamma_connmat=zeros(nnodes,nnodes,nsubs);
hnd_alpha_connmat=zeros(nnodes,nnodes,nsubs);
hnd_gamma_connmat=zeros(nnodes,nnodes,nsubs);
for ii=1:nsubs
    ii
    vis_alpha_M=zeros(nnodes,nrands);
vis_gamma_M=zeros(nnodes,nrands);
hnd_alpha_M=zeros(nnodes,nrands);
hnd_gamma_M=zeros(nnodes,nrands);

    fnamevis_full=[data_path,list_all{ii},fname_vis];
    load(fnamevis_full)
    for jj=1:size(connmat,3)%freqs
        for kk=1:size(connmat,4)%times
            aa=squeeze(connmat(:,:,jj,kk));
            aa=aa+aa'+eye(size(aa));
%             aa(zscore(aa)
% bb=tril(aa,-1);
% cc=bb;
% cc(cc>0)=zscore(bb(bb>0));
% cc(cc<1.96)=0;
            connmat(:,:,jj,kk)=aa;
        end
    end
    vis_alpha_connmat(:,:,ii)=squeeze(connmat(:,:,2,thist));    
    vis_gamma_connmat(:,:,ii)=squeeze(connmat(:,:,4,thist));
    for randcnt=1:nrands
    [vis_alpha_M(:,randcnt),q]=community_louvain(squeeze(connmat(:,:,2,thist)),gamma);
    [vis_gamma_M(:,randcnt),q]=community_louvain(squeeze(connmat(:,:,4,thist)),gamma);
    end
    D_vis_alpha = agreement(vis_alpha_M)/nrands;
    ciu_vis_alpha(:,ii,thist) = consensus_und(D_vis_alpha,tau,reps);
    D_vis_gamma = agreement(vis_gamma_M)/nrands;
    ciu_vis_gamma(:,ii,thist) = consensus_und(D_vis_gamma,tau,reps);
    
    fnamehnd_full=[data_path,list_all{ii},fname_hnd];
    load(fnamehnd_full)
    for jj=1:size(connmat,3)%freqs
        for kk=1:size(connmat,4)%times
            aa=squeeze(connmat(:,:,jj,kk));
            aa=aa+aa'+eye(size(aa));
            connmat(:,:,jj,kk)=aa;
        end
    end
    hnd_alpha_connmat(:,:,ii)=squeeze(connmat(:,:,2,thist));
%     [hnd_alpha_M(:,ii),hnd_alpha_Q(ii)]=community_louvain(squeeze(connmat(:,:,2,thist)),gamma);
    hnd_gamma_connmat(:,:,ii)=squeeze(connmat(:,:,4,thist));
%     [hnd_gamma_M(:,ii),hnd_gamma_Q(ii)]=community_louvain(squeeze(connmat(:,:,4,thist)),gamma);
    for randcnt=1:nrands
    [hnd_alpha_M(:,randcnt),q]=community_louvain(squeeze(connmat(:,:,2,thist)),gamma);
    [hnd_gamma_M(:,randcnt),q]=community_louvain(squeeze(connmat(:,:,4,thist)),gamma);
    end
    D_hnd_alpha = agreement(hnd_alpha_M)/nrands;
    ciu_hnd_alpha(:,ii,thist) = consensus_und(D_hnd_alpha,tau,reps);
    D_hnd_gamma = agreement(hnd_gamma_M)/nrands;
    ciu_hnd_gamma(:,ii,thist) = consensus_und(D_hnd_gamma,tau,reps);
end

end
% grand_sub_time_band_cond(:,:,gamma)=cat(2,reshape(ciu_vis_alpha,[nnodes,nsubs*ntimes]),reshape(ciu_vis_gamma,[nnodes,nsubs*ntimes]),reshape(ciu_hnd_alpha,[nnodes,nsubs*ntimes]),reshape(ciu_hnd_gamma,[nnodes,nsubs*ntimes]));

grand_mat_sub_time_band_cond(:,:,gcnt)=cat(2,reshape(ciu_vis_alpha,[nnodes,nsubs*ntimes]),reshape(ciu_vis_gamma,[nnodes,nsubs*ntimes]),reshape(ciu_hnd_alpha,[nnodes,nsubs*ntimes]),reshape(ciu_hnd_gamma,[nnodes,nsubs*ntimes]));
D_grand(:,:,gcnt) = agreement(squeeze(grand_mat_sub_time_band_cond(:,:,gcnt)))/squeeze(size(grand_mat_sub_time_band_cond(:,:,gcnt),2));
ciu_grand(:,gcnt) = consensus_und(squeeze(D_grand(:,:,gcnt)),tau,reps);
% thistm=2;
% grand_mat_sub_time_band_cond=squeeze(ciu_vis_alpha(:,:,thistm));%cat(2,squeeze(ciu_vis_alpha(:,:,thistm)),squeeze(ciu_vis_gamma(:,:,thistm)),squeeze(ciu_hnd_alpha(:,:,thistm)),squeeze(ciu_hnd_gamma(:,:,thistm)));
% D_grand = agreement(grand_mat_sub_time_band_cond)/size(grand_mat_sub_time_band_cond,2);
% ciu_grand = consensus_und(D_grand,tau,reps);
end
a=[];
for micnt1=1:size(ciu_grand,2)
    for micnt2=1:size(ciu_grand,2)
a(micnt1,micnt2)=mi(squeeze(ciu_grand(:,micnt1)),squeeze(ciu_grand(:,micnt2)));
    end
end
b=a./repmat(diag(a),1,size(a,2));
b=(sum(b)-1)/20;
[bm,bmw]=max(b);
[ciugs,ciugsw]=sort(ciu_grand(:,bmw));
out_file_grandmat=[out_path,'grand_mat.mat'];
grand_mat=grand_mat_sub_time_band_cond;
save(out_file_grandmat,'grand_mat')
out_file_Dgrand=[out_path,'D_grand.mat'];
save(out_file_Dgrand,'D_grand')
out_file_ciugrand=[out_path,'ciu_grand.mat'];
save(out_file_ciugrand,'ciu_grand')
wmodules=ciugsw;
finalmodules=ciugs;
out_file_modules=[out_path,'wmodules.mat'];
save(out_file_modules,'wmodules')
out_file_modules=[out_path,'finalmodules.mat'];
save(out_file_modules,'finalmodules')

