clear
clc
seeds={'lhATL','lhpIFG','lhMTL','lhAG'};%'rhATL',
bands={'Theta','Alpha','Beta','Gamma'};
% load(['/imaging/rf02/TypLexMEG/','3conn_words_',seed,'_',band,'.mat')
% load('/imaging/rf02/TypLexMEG/3conn_nonwords.mat')
% Xc=Xw; Xa=Xnw;
%Xc an Xa: 20484x2x3x17; 3 conn methods: mi, dwpli, coh#
%% over verts- bands
ntimes=2;
nsubs=17;
nverts=20484;
% which_pc=4;
nbands=length(bands);
npcs=nbands;
nseeds=length(seeds);
nconns=2;

Xc_score=zeros(nverts,ntimes,nsubs,npcs);
Xa_score=zeros(nverts,ntimes,nsubs,npcs);
Xs_score=zeros(nverts,ntimes,nsubs,npcs);

Xc_grand=zeros(nverts,ntimes,nconns,nsubs,nbands);
Xa_grand=zeros(nverts,ntimes,nconns,nsubs,nbands);
Xs_grand=zeros(nverts,ntimes,nconns,nsubs,nbands);

Xc_grand_cpc=zeros(nverts,ntimes,nsubs,nbands);
Xa_grand_cpc=zeros(nverts,ntimes,nsubs,nbands);
Xs_grand_cpc=zeros(nverts,ntimes,nsubs,nbands);

for seed_no=1:nseeds
    seed=seeds{seed_no};
    for band_no=1:nbands
        band=bands{band_no};
        load(['/imaging/rf02/TypLexMEG/','3conn_newag_orig_cncrt_',seed,'_',band,'.mat'])
        load(['/imaging/rf02/TypLexMEG/','3conn_newag_orig_abs_',seed,'_',band,'.mat'])
        Xc_grand(:,:,:,:,band_no)=squeeze(Xc(:,:,2:3,:));
        Xa_grand(:,:,:,:,band_no)=squeeze(Xa(:,:,2:3,:));
        Xs_grand(:,:,:,:,band_no)=squeeze(Xc(:,:,2:3,:))-squeeze(Xa(:,:,2:3,:));
    end
    for sub_no=1:nsubs
        sub_no
        for time_no=1:ntimes
            for band_no=1:nbands
                %concrete
                X1=squeeze(Xc_grand(:,time_no,:,sub_no,band_no));%-Xa(:,time_no,:,sub_no));
                X=X1-repmat(mean(X1),size(X1,1),1);
                [coeff,score,latent,tsquared,explained,mu]=pca(X);
                Xc_grand_cpc(:,time_no,sub_no,band_no)=score(:,1);
                %abstract
                X1=squeeze(Xa_grand(:,time_no,:,sub_no,band_no));%-Xa(:,time_no,:,sub_no));
                X=X1-repmat(mean(X1),size(X1,1),1);
                [coeff,score,latent,tsquared,explained,mu]=pca(X);
                Xa_grand_cpc(:,time_no,sub_no,band_no)=score(:,1);
                %subtract
                X1=squeeze(Xs_grand(:,time_no,:,sub_no,band_no));%-Xa(:,time_no,:,sub_no));
                X=X1-repmat(mean(X1),size(X1,1),1);
                [coeff,score,latent,tsquared,explained,mu]=pca(X);
                Xs_grand_cpc(:,time_no,sub_no,band_no)=score(:,1);
            end
        end
    end
    
    for sub_no=1:nsubs
        sub_no
        for time_no=1:ntimes
            for which_pc=1:nbands
                
                %concrete
                X1=squeeze(Xc_grand_cpc(:,time_no,sub_no,:));%-Xa(:,time_no,:,sub_no));
                X=X1-repmat(mean(X1),size(X1,1),1);
                [coeff,score,latent,tsquared,explained,mu]=pca(X);
                Xc_score(:,time_no,sub_no,:)=score;
                %abstract
                X1=squeeze(Xa_grand_cpc(:,time_no,sub_no,:));
                X=X1-repmat(mean(X1),size(X1,1),1);
                [coeff,score,latent,tsquared,explained,mu]=pca(X);
                Xa_score(:,time_no,sub_no,:)=score;
                %subtract
                X1=squeeze(Xs_grand_cpc(:,time_no,sub_no,:));
                X=X1-repmat(mean(X1),size(X1,1),1);
                [coeff,score,latent,tsquared,explained,mu]=pca(X);
                Xs_score(:,time_no,sub_no,:)=score;
            end
        end
    end
    
    out_path=['/imaging/rf02/TypLexMEG/icaanalysis_results/conn_pca/WB_PCA_cohppc_acrossverts_SL_allbands_Xcscore_',seed,'.mat'];
    save(out_path,'Xc_score')
    
    out_path=['/imaging/rf02/TypLexMEG/icaanalysis_results/conn_pca/WB_PCA_cohppc_acrossverts_SL_allbands_Xascore_',seed,'.mat'];
    save(out_path,'Xa_score')
    
    out_path=['/imaging/rf02/TypLexMEG/icaanalysis_results/conn_pca/WB_PCA_cohppc_acrossverts_SL_allbands_Xsscore_',seed,'.mat'];
    save(out_path,'Xs_score')
end

