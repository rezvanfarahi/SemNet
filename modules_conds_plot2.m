clear
clc
% close all
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
out_file_modules=[out_path,'wmodules_mi_maxctf_3conds_3bands.mat'];
load(out_file_modules)
out_file_modules=[out_path,'finalmodules_mi_maxctf_3conds_3bands.mat'];
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
            parmat(:,jj,kk)=(participation_coef(squeeze(connmat_thresh(:,:,jj,kk)),finalmodules));%zscore(participation_coef_sign(cc(wmodules,wmodules),finalmodules));%%mean(abs(squeeze(connmat_thresh(:,:,jj,kk))),2);%module_degree_zscore(squeeze(connmat_thresh(:,:,jj,kk)),finalmodules);%%gateway_coef_sign(cc(wmodules,wmodules),finalmodules,1);%participation_coef(aa2(wmodules,wmodules),finalmodules);%
        end
    end
    
    grand_connmat(ii,:,:,:,:,condcnt)=connmat1;
%     connmat_thresh(connmat_thresh>0)=1;
    grand_connmat_thresh(ii,:,:,:,:,condcnt)=connmat_thresh;
    grand_parmat(ii,:,:,:,condcnt)=parmat;
%     for fffcnt=1:size(parmat,2)
%     grand_parmat(ii,:,fffcnt,:,condcnt)=zscore(parmat(:,fffcnt,:));
%     end
    
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
                betweenness_modules(ii,:,jj,kk,ll)=strengths_und(abs(squeeze(thresh_connmat_modules(ii,:,:,jj,kk,ll))));%betweenness_wei(squeeze(thresh_connmat_modules(ii,:,:,jj,kk,ll)));
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

mean_connmat=squeeze(mean(grand_connmat_thresh,1));% thresh_connmat_modules
std_connmat=squeeze(std(thresh_connmat_modules,1));

mconnmat_alpha=squeeze(mean_connmat(:,:,2,:,:));
mconnmat_gamma=squeeze(mean_connmat(:,:,4,:,:));
sconnmat_alpha=squeeze(std_connmat(:,:,2,:,:));
sconnmat_gamma=squeeze(std_connmat(:,:,4,:,:));

thiscond=[2,3,3];
refcond=[1,1,2];
figure,
kk=0;
nfrqs=3;
centre_posname='/imaging/rf02/Semnet/BrainNetViewer_20170403/rezvan/nodeinfo_parc1_semnet.mat';
load(centre_posname)
modsig_pathin='/imaging/rf02/Semnet/stc/permutation/connectivity/parcellated/maxctf/';
modsig_fnamesin={'modules_3D_overtimeaparc_mod_Visual_Hear_Gamma_05.mat',
'modules_3D_overtimeaparc_mod_Visual_Hand_Beta_05.mat',
'modules_3D_overtimeaparc_mod_Visual_Hand_Alpha_05.mat',
'modules_3D_overtimeaparc_mod_Hear_Hand_Beta_05.mat'
};
for ii=1:nconds
    figure,
    kk=0;
    for ff=1:nfrqs
        for jj=1:ntimes
            
            kk=kk+1;
            sp1=subplot(nfrqs,ntimes,kk);
            dd=squeeze(mean_connmat(:,:,2:end,:,refcond(ii)))-squeeze(mean_connmat(:,:,2:end,:,thiscond(ii)));
            
            imagesc(squeeze(mean_connmat(:,:,ff+1,jj,refcond(ii)))-squeeze(mean_connmat(:,:,ff+1,jj,thiscond(ii))),[-0.01 0.01])
            
            colormap(hot)
            clear thist_MTpos
            if ii==1 && ff+1==4
                load([modsig_pathin,modsig_fnamesin{1}])
                thist_MTpos=squeeze(this_MTpos(:,:,jj));
            elseif ii==2 && ff+1==3
                load([modsig_pathin,modsig_fnamesin{2}])
                thist_MTpos=squeeze(this_MTpos(:,:,jj));
            elseif ii==2 && ff+1==2
                load([modsig_pathin,modsig_fnamesin{3}])
                thist_MTpos=squeeze(this_MTpos(:,:,jj));
            elseif ii==3 && ff+1==3
                load([modsig_pathin,modsig_fnamesin{4}])                
                thist_MTpos=squeeze(this_MTpos(:,:,jj));
                
            end
            if exist('thist_MTpos','var')
                ii
                jj
                ff
                thist_MTpos
%             thist_MTpos=thist_MTpos~=0;
%             for nn1=1:size(thist_MTpos,1)
%                 for nn2=1:size(thist_MTpos,2)
%                     if thist_MTpos(nn1,nn2)
%                         thist_MTpos(nn1,nn2)
%                         annotation(sp1,'rectangle',...
%                             [(nn1-1)/size(thist_MTpos,1) (nn2-1)/size(thist_MTpos,1) 1/size(thist_MTpos,1) 1/size(thist_MTpos,1)],...
%                             'Color',[0 0 0],...
%                             'LineWidth',2,...
%                             'LineStyle',':');
%                     end
%                 end
%             end
            end
            
        end
    end
end


% % % % for fn=1:length(modsig_fnamesin)
% % % %     modsig_fnamein=modsig_fnamesin{fn};
% % % % load([modsig_pathin,modsig_fnamein])%load centre_pos
% % % % % activemodules=squeeze(sum(abs(this_MTpos),2));
% % % % for tt=1:ntimes
% % % %     this_centre=centre_pos;
% % % %     switch fn
% % % %         case 1
% % % %     thismat=squeeze(grand_connmat_thresh(:,:,:,4,tt,1))-squeeze(grand_connmat_thresh(:,:,:,4,tt,2));
% % % %     strength_diff=mean(squeeze(mean_connmat(:,:,4,tt,1)),2)-mean(squeeze(mean_connmat(:,:,4,tt,2)),2);
% % % %     case 2
% % % %     thismat=squeeze(grand_connmat_thresh(:,:,:,3,tt,1))-squeeze(grand_connmat_thresh(:,:,:,3,tt,3));
% % % %     strength_diff=mean(squeeze(mean_connmat(:,:,3,tt,1)),2)-mean(squeeze(mean_connmat(:,:,3,tt,3)),2);
% % % %     case 3
% % % %     thismat=squeeze(grand_connmat_thresh(:,:,:,2,tt,1))-squeeze(grand_connmat_thresh(:,:,:,2,tt,3));
% % % %     strength_diff=mean(squeeze(mean_connmat(:,:,3,tt,1)),2)-mean(squeeze(mean_connmat(:,:,3,tt,3)),2);
% % % %     case 4
% % % %     thismat=squeeze(grand_connmat_thresh(:,:,:,3,tt,2))-squeeze(grand_connmat_thresh(:,:,:,3,tt,3));
% % % %     strength_diff=mean(squeeze(mean_connmat(:,:,3,tt,2)),2)-mean(squeeze(mean_connmat(:,:,3,tt,3)),2);
% % % %     end
% % % %     
% % % %     this_centre(:,5)=abs(strength_diff);
% % % %     cname=[out_path,modsig_fnamein(1:end-7),'_centre_pos_',timesmin{tt},'.mat'];
% % % % % % %     save(cname,'this_centre')
% % % %     for mm=1:size(this_MTpos,1)
% % % %         for nn=1:size(this_MTpos,2)
% % % %         if this_MTpos(mm,nn,tt)==0
% % % %             thismat(:,finalmodules==mm,finalmodules==nn)=0;
% % % %         end
% % % %         end
% % % %     end
% % % %     thismat_uv=zeros(nnodes,nnodes);
% % % %     for uv1=1:nnodes
% % % %         for uv2=1:nnodes
% % % %             if sum(abs(squeeze(thismat(:,uv1,uv2))))>0
% % % %             thisp=signrank(squeeze(thismat(:,uv1,uv2)));
% % % %             if thisp<=0.05
% % % %                 thismat_uv(uv1,uv2)=mean(squeeze(thismat(:,uv1,uv2)));
% % % %             end
% % % %             end
% % % %         end
% % % %     end
% % % %           thismat=thismat_uv;  
% % % %             
% % % %     thismat=thismat+thismat';
% % % %     matname=[out_path,modsig_fnamein(1:end-7),'_connmat_p005_',timesmin{tt},'.mat'];
% % % %     save(matname,'thismat')
% % % %     figure, imagesc(thismat)
% % % % end
% % % % end


for fn=1:length(modsig_fnamesin)
    modsig_fnamein=modsig_fnamesin{fn};
load([modsig_pathin,modsig_fnamein])%load centre_pos
% activemodules=squeeze(sum(abs(this_MTpos),2));
for tt=1:ntimes
    this_centre=centre_pos;
    switch fn
        case 1
    thismat=squeeze(mean_connmat(:,:,4,tt,1))-squeeze(mean_connmat(:,:,4,tt,2));
    strength_diff=mean(squeeze(mean_connmat(:,:,4,tt,1)),2)-mean(squeeze(mean_connmat(:,:,4,tt,2)),2);
    case 2
    thismat=squeeze(mean_connmat(:,:,3,tt,1))-squeeze(mean_connmat(:,:,3,tt,3));
    strength_diff=mean(squeeze(mean_connmat(:,:,3,tt,1)),2)-mean(squeeze(mean_connmat(:,:,3,tt,3)),2);
    case 3
    thismat=squeeze(mean_connmat(:,:,2,tt,1))-squeeze(mean_connmat(:,:,2,tt,3));
    strength_diff=mean(squeeze(mean_connmat(:,:,3,tt,1)),2)-mean(squeeze(mean_connmat(:,:,3,tt,3)),2);
    case 4
    thismat=squeeze(mean_connmat(:,:,3,tt,2))-squeeze(mean_connmat(:,:,3,tt,3));
    strength_diff=mean(squeeze(mean_connmat(:,:,3,tt,2)),2)-mean(squeeze(mean_connmat(:,:,3,tt,3)),2);
    end
    this_centre(:,5)=abs(strength_diff);
    
    for mm=1:size(this_MTpos,1)
        for nn=1:size(this_MTpos,2)
        if this_MTpos(mm,nn,tt)==0
            thismat(finalmodules==mm,finalmodules==nn)=0;
        end
        end
    end
    thisabs=abs(thismat);
    thismat(thisabs<max(thisabs(:))/2)=0;    
    thismat=thismat+thismat';
%     strength_diff=mean(thisabs+thisabs',2);
%     if sum(strength_diff)>0
%     this_centre(:,5)=abs(strength_diff+min(strength_diff(strength_diff>0)));
%     else
%         this_centre(:,5)=1;
%     end
    cname=[out_path,modsig_fnamein(1:end-7),'_centre_pos_',timesmin{tt},'.mat'];
% % %     save(cname,'this_centre')
    matname=[out_path,modsig_fnamein(1:end-7),'_connmat_halfmax_',timesmin{tt},'.mat'];
%     save(matname,'thismat')
%     figure, imagesc(thismat)
end
end

% colormap(jet)
% colorbar

% figure,
% kk=0;
% for ii=1:nconds
%     for jj=1:ntimes
%         kk=kk+1;
%         subplot(ntimes,nconds,kk)
%         bar(squeeze(sum(mconnmat_alpha(:,:,jj,refcond(ii)),2))-squeeze(sum(mconnmat_alpha(:,:,jj,thiscond(ii)),2)))%,[-0.01 0.01])
%     end
% end




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





% figure,
% kk=0;
% for ii=1:nconds
%     for jj=1:ntimes
%         kk=kk+1;
%         subplot(ntimes,nconds,kk)
%         imagesc(squeeze(mconnmat_gamma(:,:,jj,refcond(ii)))-squeeze(mconnmat_gamma(:,:,jj,thiscond(ii))),[-0.01 0.01])
%     end
% end
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

% figure,
% kk=0;
% for ii=1:nconds
%     for jj=1:ntimes
%         kk=kk+1;
%         subplot(ntimes,nconds,kk)
%         bar(squeeze(sum(mconnmat_gamma(:,:,jj,refcond(ii)),2))-squeeze(sum(mconnmat_gamma(:,:,jj,thiscond(ii)),2)))%,[-0.01 0.01])
%         
%     end
% end

% figure,
% kk=0;
% for ii=1:nconds
%     for jj=1:ntimes
%         kk=kk+1;
%         subplot(ntimes,nconds,kk)
%         imagesc(squeeze(mconnmat_alpha(:,:,jj,ii)))
%     end
% end
% colormap(hot)

% figure,
% kk=0;
% nfrqs=3;
% for ii=1:nconds
%     figure,
%     kk=0;
%     for ff=1:nfrqs
%         for jj=1:ntimes
%             
%             kk=kk+1;
%             subplot(nfrqs,ntimes,kk)
%             imagesc(squeeze(mean_connmat(:,:,ff+1,jj,ii)));%imagesc(squeeze(mean_connmat(:,:,ff+1,jj,refcond(ii)))-squeeze(mean_connmat(:,:,ff+1,jj,thiscond(ii))),[-0.01 0.01])
%         end
%     end
%     colormap(hot)
% 
% end
% mn=squeeze(mean(mean_connmat,5));
% mn=squeeze(mean(mn,4));
% mn=squeeze(mean(mn(:,:,2:end),3));
% figure, 
% createfigure(mn)
        
% figure,
% kk=0;
% for ii=1:nconds
%     for jj=1:ntimes
%         kk=kk+1;
%         subplot(ntimes,nconds,kk)
%         imagesc(squeeze(mconnmat_gamma(:,:,jj,ii)))
%     end
% end       
% colormap(hot)
        
load('/imaging/rf02/Semnet/graph_analysis/labels_aparc_mod.mat')
lmod=label_names(wmodules,:);
% % % % out_file_grandconnmat=[out_path,'grand_connmat_mi_maxctf_modules3conds_3bands.mat'];
% % % % save(out_file_grandconnmat,'grand_connmat')
% % % % out_file_grandconnmat=[out_path,'grand_connmat_mi_maxctf_thresh_modules3conds_3bands.mat'];
% % % % save(out_file_grandconnmat,'grand_connmat_thresh')
% % % % out_file_grandconnmat=[out_path,'grand_participation_mi_maxctf_nodes_modules3conds_3bands.mat'];
% % % % save(out_file_grandconnmat,'grand_parmat')
% % % % out_file_grandconnmat=[out_path,'grand_modconnmat_mi_maxctf_thresh_modules3conds_3bands.mat'];
% % % % save(out_file_grandconnmat,'thresh_connmat_modules')

 
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
pval_par=ones(nnodes,ntimes,nbands);
thisev=1;
refev=3;
for thisfreq=1:4;
pval_mat_par=zeros(nnodes,ntimes);
for ii=1:size(pval_mat_par,1)
    for jj=1:size(pval_mat_par,2)
            pval_mat_par(ii,jj)=signrank(squeeze(grand_parmat(:,ii,thisfreq,jj,refev))-squeeze(grand_parmat(:,ii,thisfreq,jj,thisev)));
    end
end
pval_par1=pval_mat_par;
pval_par1(pval_par1>0.005)=1;
pval_par(:,:,thisfreq)=pval_par1;
for ii=1:ntimes
FDR(:,ii,thisfreq) = mafdr(pval_mat_par(:,ii));
end
FDR(FDR>=0.05)=1;
min(FDR(:))
end
% pval_mat_par(pval_mat_par>0.05)=1;
            
close all            

pthresh=0.025;
thisev=3;
refev=2;
nperm=1000;
freqs=3;
nbands=length(freqs);
tval_par=zeros(nnodes,ntimes,nbands);
pval_par=ones(nnodes,ntimes,nbands);
fcn=0;
% pt=mean(grand_parmat(:))+std(grand_parmat(:));
% grand_parmat1=grand_parmat(:,:,2:end,:,:);
% grand_parmat1(grand_parmat1<pt)=0;
% grand_parmat1(grand_parmat1>=pt)=1;
% grand_parmat(:,:,2:end,:,:)=grand_parmat1;

for thisfreq=freqs;
    fcn=fcn+1;
pval_mat_par=zeros(nnodes,ntimes);
for ii=1:size(pval_mat_par,1)
    for jj=1:size(pval_mat_par,2)
        thissubt=squeeze(grand_parmat(:,ii,thisfreq,jj,refev))-squeeze(grand_parmat(:,ii,thisfreq,jj,thisev));
            [pval_mat_par(ii,jj),h,stat1]=signrank(thissubt);%[h,pval_mat_par(ii,jj),ci,stat1]=ttest(thissubt);
            tval_par(ii,jj,fcn)=stat1.zval;
    end
end

pval_par(:,:,fcn)=pval_mat_par;

end            
tval_par(pval_par>=pthresh)=0;
tval_orig=tval_par;
tval_sum=squeeze(sum(sum(abs(tval_par),3),2));

tval_max_perm=zeros(nperm,1);
for np=1:nperm
    np
    r=sign(rand(nsubs,1)-0.5);
    tval_par=zeros(nnodes,ntimes,nbands);
    pval_par=ones(nnodes,ntimes,nbands);
    fcn=0;
    for thisfreq=freqs
        fcn=fcn+1;
        pval_mat_par=zeros(nnodes,ntimes);
        for ii=1:size(pval_mat_par,1)
            for jj=1:size(pval_mat_par,2)
                thissubt=squeeze(grand_parmat(:,ii,thisfreq,jj,refev))-squeeze(grand_parmat(:,ii,thisfreq,jj,thisev));
                [pval_mat_par(ii,jj),h,stat1]=signrank(thissubt.*r);
                tval_par(ii,jj,fcn)=stat1.zval;
            end
        end
        
        pval_par(:,:,fcn)=pval_mat_par;
        
    end
    tval_par(pval_par>=pthresh)=0;
    tval_sum_perm=squeeze(sum(sum(abs(tval_par),3),2));
    tval_max_perm(np)=max(tval_sum_perm);
end
pval_perm=zeros(size(tval_sum));
for pn=1:length(pval_perm)
    pval_perm(pn)=1-(sum(tval_sum(pn)>tval_max_perm)/nperm);
end
lmod(pval_perm<0.05)
    
    
    
    
    
    


