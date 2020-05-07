clear
clc
addpath('/imaging/rf02/Semnet/NBS1.2/')
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
vis_alpha_connmat=zeros(74,74,nsubs);
vis_gamma_connmat=zeros(74,74,nsubs);
hnd_alpha_connmat=zeros(74,74,nsubs);
hnd_gamma_connmat=zeros(74,74,nsubs);

% hnd_connmat=zeros(74,74,2,nsubs);
timesmin={'50ms','150ms','250ms'};
timesmax={'250ms','350ms','450ms'};
thist=2;
for ii=1:nsubs
    ii
    fnamevis_full=[data_path,list_all{ii},fname_vis];
    load(fnamevis_full)
    for jj=1:size(connmat,3)%freqs
        for kk=1:size(connmat,4)%times
            aa=squeeze(connmat(:,:,jj,kk));
            aa=aa+aa'+eye(size(aa));
            connmat(:,:,jj,kk)=aa;
        end
    end
    vis_alpha_connmat(:,:,ii)=squeeze(connmat(:,:,2,thist));
    vis_gamma_connmat(:,:,ii)=squeeze(connmat(:,:,4,thist));
    
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
    hnd_gamma_connmat(:,:,ii)=squeeze(connmat(:,:,4,thist));
end
grand_connmat=cat(3,vis_alpha_connmat,vis_gamma_connmat,hnd_alpha_connmat,hnd_gamma_connmat);
fname_out=[data_path,'/NBS1.2/Rezvan/','Matrix_vishnd_alphagamma_',timesmin{thist},'_',timesmax{thist},'.mat'];
save(fname_out,'grand_connmat')
gconnmat2=cat(3,vis_gamma_connmat,hnd_gamma_connmat);
fname_out=[data_path,'/NBS1.2/Rezvan/','Matrix_vishnd_gamma_',timesmin{thist},'_',timesmax{thist},'.mat'];
save(fname_out,'gconnmat2')

gconnmat3=cat(3,vis_alpha_connmat,hnd_alpha_connmat);%vis_alpha_connmat-hnd_alpha_connmat;%
fname_out=[data_path,'/NBS1.2/Rezvan/','Matrix_vishnd_alpha_',timesmin{thist},'_',timesmax{thist},'.mat'];
save(fname_out,'gconnmat3')


dmat=cat(1,eye(nsubs),eye(nsubs));
dmat=[dmat,[ones(nsubs,1);-1*ones(nsubs,1)]];
% dmat=ones(nsubs,1);
fname_out=[data_path,'/NBS1.2/Rezvan/','designMatrix_vishnd_gamma_',timesmin{thist},'_',timesmax{thist},'.mat'];
save(fname_out,'dmat')
load('/imaging/rf02/Semnet/NBS1.2/Rezvan/nodeinfo_parc1.mat')
nodenames={};
for nc=1:size(centre_pos,1)
    if mod(nc,2)==1 %left
    bb=['L_',num2str(((nc+1)/2))];
    else
    bb=['R_',num2str((nc/2))];
    end
    nodenames{nc,1}=bb;
end

cont=vis_alpha_connmat-hnd_alpha_connmat;
nnodes=size(centre_pos,1);
pval=zeros(nnodes,nnodes);
for p1=1:nnodes
    for p2=1:nnodes
       [h, pval(p1,p2)]=ttest(squeeze(cont(p1,p2,:)));
    end
end


    
    
    
    

    
    
    
    
    
    
    
    
    
    