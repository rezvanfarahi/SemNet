mcc=0;
mc_gmat=zeros(36,5,12);
for nseeds=[3,5,10,15]
    nseeds
    seedcnt=seedcnt+1;
    concnt=0;
    for nconns=[0.25,0.5,1]
        mcc=mcc+1;
        nconns
        concnt=concnt+1;
        for which_parc=[1,2,3,4,5]

            nsubs=17;
            which_parc
            % nseeds=5;
            % nconns=1;
            %both SNRs
            switch which_parc
                case 1
                    nnodes=68;
                    inpath=['/imaging/rf02/parcellation/extended_simulations/funcseeds_avgsims70_offset0_seeds',num2str(nseeds),'_conns',num2str(nconns),'_bothsnrs_extglasser/aparc/'];
                    cd(inpath)
                    load ('perm_pval_parc0.mat');
%                     load ('ideal_pval_parc0.mat');
                case 2
                    nnodes=74;
                    inpath=['/imaging/rf02/parcellation/extended_simulations/funcseeds_avgsims70_offset0_seeds',num2str(nseeds),'_conns',num2str(nconns),'_bothsnrs_extglasser/aparc_mod/'];
                    cd(inpath)
                    load ('perm_pval_parc1.mat');
%                     load ('ideal_pval_parc1.mat');
                case 3
                    nnodes=148;
                    inpath=['/imaging/rf02/parcellation/extended_simulations/funcseeds_avgsims70_offset0_seeds',num2str(nseeds),'_conns',num2str(nconns),'_bothsnrs_extglasser/aparc2009/'];
                    cd(inpath)
                    load ('perm_pval_parc2.mat');
%                     load ('ideal_pval_parc2.mat');
                case 4
                    nnodes=74;
                    inpath=['/imaging/rf02/parcellation/extended_simulations/funcseeds_avgsims70_offset0_seeds',num2str(nseeds),'_conns',num2str(nconns),'_bothsnrs_extglasser/aparc2009_mod/'];
                    cd(inpath)
                    load ('perm_pval_parc3.mat');
%                     load ('ideal_pval_parc3.mat');
                case 5
                    nnodes=70;
                    inpath=['/imaging/rf02/parcellation/extended_simulations/funcseeds_avgsims70_offset0_seeds',num2str(nseeds),'_conns',num2str(nconns),'_bothsnrs_extglasser/RG/'];
                    cd(inpath)
                    load ('perm_pval_parc4.mat');
%                     load ('ideal_pval_parc4.mat');
            end
                        
            load(['/imaging/rf02/parcellation/extended_simulations/funcseeds_avgsims70_offset0_seeds',num2str(nseeds),'_conns',num2str(nconns),'_bothsnrs_extglasser/allmissedcons.mat'])
%             missedcons=squeeze(mean(missedcons,3));
            mc_gmat(:,:,mcc)=missedcons;
%             missedcons=missedcons(:,which_parc);
        end
    end
end
