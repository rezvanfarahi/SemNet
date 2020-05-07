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
fname_conds={'pacel_connmat_maxctf_Visual_aparc_mod.mat','pacel_connmat_maxctf_Hand_aparc_mod.mat','pacel_connmat_maxctf_Hear_aparc_mod.mat'};
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
out_file_modules=[out_path,'wmodules_maxctf_3conds.mat'];
load(out_file_modules)
out_file_modules=[out_path,'finalmodules_maxctf_3conds.mat'];
load(out_file_modules)
connmat1=zeros(nnodes,nnodes,nbands,ntimes);
node_perc=[0.05:0.05:0.75];

Deg_b=zeros(nsubs,nnodes,length(node_perc));

tr=zeros(nsubs,length(node_perc));
trb=zeros(nsubs,length(node_perc));

for cnt=1:nsubs
%     cnt

for condcnt=1:nconds
    fnamecond_full=[data_path,list_all{cnt},fname_conds{condcnt}];
    load(fnamecond_full)
    
    for jj=1:size(connmat,3)%freqs
        for kk=1:size(connmat,4)%times
            aa=squeeze(connmat(:,:,jj,kk));
            aa=aa+aa';%+eye(size(aa));
%             aa(zscore(aa)
% bb=tril(aa,-1);
% cc=bb;
% cc(cc>0)=zscore(bb(bb>0));
% cc(cc<1.96)=0;
            connmat1(:,:,jj,kk)=aa(wmodules,wmodules);
        end
    end
    grand_connmat(cnt,:,:,:,:,condcnt)=connmat1;
    
end
plv_matb=squeeze(mean(grand_connmat(cnt,:,:,:,:,:),6));
thresh_b=zeros(size(plv_matb));
    deg_b=zeros(size(plv_matb,2),size(plv_matb,3),size(plv_matb,4));
    
      
    plrb=plv_matb(:);
    [aab,bbb]=sort(plrb,'descend'); 
    aab2=aab(aab>0);
    
    cofi=0;
    for Cof=node_perc
        Cof
        cofi=cofi+1;
%         tr(cnt,cofi)=aa(round(Cof*length(aa)));
        trb(cnt,cofi)=aab2(round(Cof*length(aab2)));
    
        for mm=1:size(plv_matb,3)
            mm
            %         pla=squeeze(abs(plv_mata(:,:,mm,:)));
            %         tra=mean(pla(pla~=0));
            %
            %         plc=squeeze(abs(plv_matc(:,:,mm,:)));
            %         plr=[reshape(plc,1,size(plc,1)*size(plc,2)*size(plc,3)),reshape(pla,1,size(pla,1)*size(pla,2)*size(pla,3))];
            %         [aa,bb]=sort(plr,'descend');
            %         tr=aa(round(length(aa)/4));
            %         tr
            %         trc=mean(plc(plc~=0));
            %tr=mean([tra,trc]);
            for ll=1:size(plv_matb,4)                
                thresh_b(:,:,mm,ll)=threshold_absolute(plv_matb(:,:,mm,ll),trb(cnt,cofi));
                deg_b(:,cofi,mm,ll)=degrees_und(thresh_b(:,:,mm,ll));
                
                %             thresh_s(:,:,mm,ll)=threshold_absolute(plv_mats(:,:,mm,ll),0.5);
                %             deg_s(:,mm,ll)=degrees_und(thresh_s(:,:,mm,ll));
            end
        end
% %         deg_cr=zeros(size(deg_c,1),size(deg_c,2),size(deg_c,3),length([6:5:ll-5]));
% %         deg_ar=zeros(size(deg_a,1),size(deg_a,2),size(deg_a,3),length([6:5:ll-5]));
% %         cc=0;
% %         for bc=6:5:ll-5
% %             cc=cc+1;
% %             b1=bc-5; b2=b1+10;
% %             deg_cr(:,:,:,cc)=mean(deg_c(:,:,:,b1:b2),4);
% %             deg_ar(:,:,:,cc)=mean(deg_a(:,:,:,b1:b2),4);
% %         end
        
        deg_br=deg_b;
        
        deg_b1=mean(deg_br,4);deg_b2=mean(deg_b1,3);
        Deg_b(cnt,:,cofi)=deg_b2(:,cofi);
    end
end


Deg_m=Deg_b;
Degarat=Deg_m./repmat(mean(Deg_m,2),1,size(Deg_m,2),1);
figure,imagesc(node_perc*100,1:nnodes,squeeze(mean(Degarat,1))); colorbar%,ax1=subplot(1,2,1);colormap(ax1,jet),
degatest=abs(diff(Degarat,1,3))-0.05;
tmap=zeros(size(degatest,2),size(degatest,3));
for ii=1:size(tmap,1)
for jj=1:size(tmap,2)
tmap(ii,jj)=ttest(degatest(:,ii,jj));
end
end
figure,imagesc(node_perc(2:end)*100,1:nnodes,tmap), colormap(gray)%ax2=subplot(1,2,2);