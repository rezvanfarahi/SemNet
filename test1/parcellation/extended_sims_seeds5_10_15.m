clear
clc
close all
%%aparc
%SNR3-chunck1
TP_grand=zeros(4,3,5,9);%the 9: COH:leak3,noise3,leak1,noise1, imCOH:leak3,noise3,ideal,leak1,noise1,
FP_grand=zeros(4,3,5,9);
TN_grand=zeros(4,3,5,9);
FN_grand=zeros(4,3,5,9);
seedcnt=0;
for nseeds=[10]
    nseeds
    seedcnt=seedcnt+1;
    concnt=0;
    for nconns=[1]
        nconns
        concnt=concnt+1;
        for which_parc=[1,2,3,4,5]
            nsims=35;%length(FileName);
            nsubs=17;
            which_parc
            % nseeds=5;
            % nconns=1;
            %both SNRs
            switch which_parc
                case 1
                    nnodes=68;
                    inpath=['/imaging/rf02/parcellation/extended_simulations/sims70_offset0_seeds',num2str(nseeds),'_conns',num2str(nconns),'_bothsnrs_ext25/aparc/'];
                    cd(inpath)
                    aa=dir ('*parc_connmat*.mat');
                case 2
                    nnodes=74;
                    inpath=['/imaging/rf02/parcellation/extended_simulations/sims70_offset0_seeds',num2str(nseeds),'_conns',num2str(nconns),'_bothsnrs_ext25/aparc_mod/'];
                    cd(inpath)
                    aa=dir ('*parc_connmat*.mat');
                case 3
                    nnodes=148;
                    inpath=['/imaging/rf02/parcellation/extended_simulations/sims70_offset0_seeds',num2str(nseeds),'_conns',num2str(nconns),'_bothsnrs_ext25/aparc2009/'];
                    cd(inpath)
                    aa=dir ('*parc_connmat*.mat');
                case 4
                    nnodes=74;
                    inpath=['/imaging/rf02/parcellation/extended_simulations/sims70_offset0_seeds',num2str(nseeds),'_conns',num2str(nconns),'_bothsnrs_ext25/aparc2009_mod/'];
                    cd(inpath)
                    aa=dir ('*parc_connmat*.mat');
                case 5
                    nnodes=70;
                    inpath=['/imaging/rf02/parcellation/extended_simulations/sims70_offset0_seeds',num2str(nseeds),'_conns',num2str(nconns),'_bothsnrs_ext25/RG/'];
                    cd(inpath)
                    aa=dir ('*parc_connmat*.mat');
            end
            
            
            load(['/imaging/rf02/parcellation/extended_simulations/sims70_offset0_seeds',num2str(nseeds),'_conns',num2str(nconns),'_bothsnrs_ext25/missedcons/allmissedcons.mat'])
            missedcons=squeeze(mean(missedcons,3));
            missedcons=missedcons(:,which_parc);
            if nseeds==3 && nconns==1
                nseeds%load('connmat_35sims.mat')
            else
%                 [FileName,PathName,FilterIndex] = uigetfile('*.mat',['Select parc',num2str(which_parc),'seed',num2str(nseeds),'conn',num2str(nconns)],'MultiSelect','on');
                connmat=zeros(nnodes,nnodes,10,nsubs,nsims);%leakageSNR3,noiseSNR3,ideal,leakageSNR1,noiseSNR1 (COH 1st, imCOH2nd)
                cc=0;
                for fcnt=1:35%nsims
                    cc=cc+1;
                    a=load([inpath,aa(fcnt).name]);%load([PathName,FileName{cc}]);
                    switch which_parc
                        case 1
                    connmat(:,:,:,:,fcnt)=a.parc0_connmat;%(:,:,:,:,fcnt)
                    clear a
                        case 2
                            connmat(:,:,:,:,fcnt)=a.parc1_connmat;
                            clear a
                        case 3
                            connmat(:,:,:,:,fcnt)=a.parc2_connmat;
                            clear a
                        case 4
                            connmat(:,:,:,:,fcnt)=a.parc3_connmat;
                            clear a
                        case 5
                            connmat(:,:,:,:,fcnt)=a.parc4_connmat;
                            clear a
                    end
                end
                
                %null
                switch which_parc
                    case 1
                        inpath='/imaging/rf02/parcellation/extended_simulations/null_sims40_offset0_bothsnrs/aparc/';
                        cd(inpath)
                        aa=dir ('*parc_connmat*.mat');
                    case 2
                        inpath='/imaging/rf02/parcellation/extended_simulations/null_sims40_offset0_bothsnrs/aparc_mod/';
                        cd(inpath)
                        aa=dir ('*parc_connmat*.mat');
                    case 3
                        inpath='/imaging/rf02/parcellation/extended_simulations/null_sims40_offset0_bothsnrs/aparc2009/';
                        cd(inpath)
                        aa=dir ('*parc_connmat*.mat');
                    case 4
                        inpath='/imaging/rf02/parcellation/extended_simulations/null_sims40_offset0_bothsnrs/aparc2009_mod/';
                        cd(inpath)
                        aa=dir ('*parc_connmat*.mat');
                    case 5
                        inpath='/imaging/rf02/parcellation/extended_simulations/null_sims40_offset0_bothsnrs/RG/';
                        cd(inpath)
                        aa=dir ('*parc_connmat*.mat');
                end
                
                nullconnmat=zeros(nnodes,nnodes,8,nsubs,nsims);
%                 [FileName,PathName,FilterIndex] = uigetfile('*.mat',['Select parc',num2str(which_parc),'null'],'MultiSelect','on');
                cc=0;
                for fcnt=1:35%nsims
                    cc=cc+1;
                    a=load([inpath,aa(fcnt).name]);%load([PathName,FileName{cc}]);
                    switch which_parc
                        case 1
                            nullconnmat(:,:,:,:,fcnt)=a.parc0_connmat;%(:,:,:,:,fcnt)
                            clear a
                        case 2
                            nullconnmat(:,:,:,:,fcnt)=a.parc1_connmat;
                            clear a
                        case 3
                            nullconnmat(:,:,:,:,fcnt)=a.parc2_connmat;
                            clear a
                        case 4
                            nullconnmat(:,:,:,:,fcnt)=a.parc3_connmat;
                            clear a
                        case 5
                            nullconnmat(:,:,:,:,fcnt)=a.parc4_connmat;
                            clear a
                    end
                end
                nullconnmat=mean(abs(nullconnmat),5);
                base_connmat=repmat(nullconnmat(:,:,[2,1,1,4,3,6,5,5,8,7],:),[1,1,1,1,nsims]);
                clear nullmat
                % pre_base_connmat=base_connmat;
                % pre_connmat=connmat;
                % for cc1=1:10
                %     for cc2=1:nsubs
                %     for cc3=1:nsims
                %
                %         thisbase=abs(base_connmat(:,:,cc1,cc2,cc3)+triu(1000*ones(nnodes,nnodes)));
                %         thisconn=abs(connmat(:,:,cc1,cc2,cc3)+triu(1000*ones(nnodes,nnodes)));
                %
                %         thisbase=(thisbase-min(thisbase(thisbase<1000)))/(max(thisbase(thisbase<1000))-min(thisbase(thisbase<1000)));
                %         thisconn=(thisconn-min(thisconn(thisconn<1000)))/(max(thisconn(thisconn<1000))-min(thisconn(thisconn<1000)));
                %         connmat(:,:,cc1,cc2,cc3)=tril(thisconn,-1);
                %         base_connmat(:,:,cc1,cc2,cc3)=tril(thisbase,-1);
                %     end
                % end
                % end
                % min(connmat(:))
                basecorr_connmat=abs(connmat)-abs(base_connmat);%ideal 3COH,8imCOH
                pval=zeros(nnodes,nnodes,10,nsims);
                disp('ttest')
                for cc1=1:nnodes
%                     cc1
                    for cc2=1:nnodes
%                         cc2
                        if cc1>cc2
                            for cc3=1:10
                                for cc4=1:nsims
                                    thisvect=squeeze(basecorr_connmat(cc1,cc2,cc3,:,cc4));
                                    [h,pval(cc1,cc2,cc3,cc4),ci]=ttest(thisvect,0,[],1);
                                end
                            end
                        end
                    end
                end
                pfdrmat=zeros(size(pval));
                for cc1=1:10
                    cc1
                    for cc2=1:nsims
                        cc2
                        pp=pval(:,:,cc1,cc2)+triu(1000*ones(nnodes,nnodes));
                        pp=pp(pp<1000);
                        pfdr=mafdr(pp);
                        ii=0;
                        for cc3=1:nnodes
                            for cc4=1:nnodes
                                if cc3<cc4
                                    ii=ii+1;
                                    pfdrmat(cc4,cc3,cc1,cc2)=pfdr(ii);
                                end
                            end
                        end
                    end
                end
                ntests=nnodes*(nnodes-1)/2;
                % pthresh=pval/ntests;
                disp('significant conns')
                
                %%% finding threshold
% %                 flag=ntests-1;
% %                 flcnt=-1;
% %                 while thisflag>0
% %                     flcnt=flcnt+1;
% %                 pfdrmat=pval*(ntests-flcnt);
% %                 for cc=1:nsims
% %                     thismat=squeeze(pfdrmat(:,:,3,cc));
% %                         thismat(thismat>=0.05)=-1;
% %                         thismat(thismat>-1)=0;
% %                       
% %                             thismat=thismat+1;
% %                         
% %                         thismat=tril(thismat,-1);
% %                         if sum(thismat(:))~=ceil(nconns*nseeds*(nseeds-1)/2)-missedcons(cc)
% %                         pfdrmat_final(:,:,cc1,cc2)=thismat;
                    
                pfdrmat_final=zeros(size(pfdrmat));
                bcm=squeeze(mean(basecorr_connmat(:,:,3,:,:),4));
                for cc1=1:10
%                     cc1
                    for cc2=1:nsims
%                         cc2
                        thismat=squeeze(pfdrmat(:,:,cc1,cc2));
                        qq=bcm(:,:,cc2)+triu(1000*ones(nnodes,nnodes));
                        thismat(thismat>=0.05)=-1;
                        this
                        thismat(thismat>-1)=0;
                        if cc1==3 %|| cc1==8 %ideal cases
                            thismat=thismat+1;
                        else
                            thismat=(thismat+2)*2;
                        end
                        thismat=tril(thismat,-1);
                        pfdrmat_final(:,:,cc1,cc2)=thismat;
                    end
                end
                eval_mat=pfdrmat_final(:,:,[1,2,4,5,6,7,8,9,10],:)- pfdrmat_final(:,:,[3,3,3,3,3,3,3,3,3],:);
                TP=zeros(9,nsims);
                TN=zeros(9,nsims);
                FP=zeros(9,nsims);
                FN=zeros(9,nsims);
                for cc1=1:9
                    for cc2=1:nsims
                        thiseval=eval_mat(:,:,cc1,cc2);
                        refvect=pfdrmat_final(:,:,3,cc2)+triu(1000*ones(nnodes,nnodes));
                        msdcn=missedcons(cc2);
                        refpos=max([sum(sum(refvect==1)),1])+msdcn;
                        refneg=max([sum(sum(refvect==0)),1]);
                        TP(cc1,cc2)=sum(sum(thiseval==3))/refpos;
                        FP(cc1,cc2)=sum(sum(thiseval==4))/refpos;
                        TN(cc1,cc2)=sum(sum(thiseval==2))/refneg;
                        FN(cc1,cc2)=(sum(sum(thiseval==1))+msdcn)/refpos;
                    end
                end
                
                TP_grand(seedcnt,concnt,which_parc,:)=mean(TP,2); %the 9: COH:leak3,noise3,leak1,noise1, imCOH:leak3,noise3,ideal,leak1,noise1,
                FP_grand(seedcnt,concnt,which_parc,:)=mean(FP,2);
                TN_grand(seedcnt,concnt,which_parc,:)=mean(TN,2);
                FN_grand(seedcnt,concnt,which_parc,:)=mean(FN,2);
                
            end
        end
        output_path='/imaging/rf02/parcellation/extended_simulations/TP_grand.mat';
        save(output_path,'TP_grand')
        output_path='/imaging/rf02/parcellation/extended_simulations/TN_grand.mat';
        save(output_path,'TN_grand')
        output_path='/imaging/rf02/parcellation/extended_simulations/FP_grand.mat';
        save(output_path,'FP_grand')
        output_path='/imaging/rf02/parcellation/extended_simulations/FN_grand.mat';
        save(output_path,'FN_grand')
    end
end


% pfdrmat(pfdrmat>=0.05)=-13
% pfdrmat(pfdrmat>-1)=0;
% pfdrmat0=pfdrmat+1;
%

% save('connmat_35sims.mat','connmat')
% %
% % %%aparc_mod
% % %SNR3-chunck1
% % nsims=35;%length(FileName);
% % nnodes=74;
% % nsubs=17;
% % connmat=zeros(nnodes,nnodes,10,nsubs,nsims);%leakageSNR3,noiseSNR3,ideal,leakageSNR1,noiseSNR1 (COH 1st, imCOH2nd)
% % %both SNRs
% % cd(['/imaging/rf02/parcellation/extended_simulations/sims70_offset0_seeds',num2str(nseeds),'_conns',num2str(nconns),'_bothsnrs_ext25/aparc_mod'])
% % if nseeds==3 && nconns==1
% %     load('connmat_35sims.mat')
% % else
% % [FileName,PathName,FilterIndex] = uigetfile('*.mat','Select the Mat-file','MultiSelect','on');
% % cc=0;
% % for fcnt=1:35%nsims
% %     cc=cc+1;
% % a=load([PathName,FileName{cc}]);
% % connmat(:,:,:,:,fcnt)=a.parc1_connmat;%(:,:,:,:,fcnt)
% % end
% % end
% % % save('connmat_35sims.mat','connmat')
% %
% % %%aparc2009
% % %SNR3-chunck1
% % nsims=35;%length(FileName);
% % nnodes=148;
% % nsubs=17;
% % connmat=zeros(nnodes,nnodes,10,nsubs,nsims);%leakageSNR3,noiseSNR3,ideal,leakageSNR1,noiseSNR1 (COH 1st, imCOH2nd)
% % %both SNRs
% % cd(['/imaging/rf02/parcellation/extended_simulations/sims70_offset0_seeds',num2str(nseeds),'_conns',num2str(nconns),'_bothsnrs_ext25/aparc2009'])
% % if nseeds==3 && nconns==1
% %     load('connmat_35sims.mat')
% % else
% % [FileName,PathName,FilterIndex] = uigetfile('*.mat','Select the Mat-file','MultiSelect','on');
% % cc=0;
% % for fcnt=1:35%nsims
% %     cc=cc+1;
% % a=load([PathName,FileName{cc}]);
% % connmat(:,:,:,:,fcnt)=a.parc2_connmat;%(:,:,:,:,fcnt)
% % end
% % end
% % % save('connmat_35sims.mat','connmat')
% %
% % %%aparc2009_mod
% % %SNR3-chunck1
% % nsims=35;%length(FileName);
% % nnodes=74;
% % nsubs=17;
% % connmat=zeros(nnodes,nnodes,10,nsubs,nsims);%leakageSNR3,noiseSNR3,ideal,leakageSNR1,noiseSNR1 (COH 1st, imCOH2nd)
% %
% % %both SNRs
% % cd(['/imaging/rf02/parcellation/extended_simulations/sims70_offset0_seeds',num2str(nseeds),'_conns',num2str(nconns),'_bothsnrs_ext25/aparc2009_mod'])
% % if nseeds==3 && nconns==1
% %     load('connmat_35sims.mat')
% % else
% % [FileName,PathName,FilterIndex] = uigetfile('*.mat','Select the Mat-file','MultiSelect','on');
% % cc=0;
% % for fcnt=1:35%nsims
% %     cc=cc+1;
% % a=load([PathName,FileName{cc}]);
% % connmat(:,:,:,:,fcnt)=a.parc3_connmat;%(:,:,:,:,fcnt)
% % end
% % end
% % % save('connmat_35sims.mat','connmat')
% %
% % %%RG
% %
% % nsims=35;%length(FileName);
% % nnodes=70;
% % nsubs=17;
% % connmat=zeros(nnodes,nnodes,10,nsubs,nsims);%leakageSNR3,noiseSNR3,ideal,leakageSNR1,noiseSNR1 (COH 1st, imCOH2nd)
% %
% % %both SNRs
% % cd(['/imaging/rf02/parcellation/extended_simulations/sims70_offset0_seeds',num2str(nseeds),'_conns',num2str(nconns),'_bothsnrs_ext25/RG'])
% % if nseeds==3 && nconns==1
% %     load('connmat_35sims.mat')
% % else
% % [FileName,PathName,FilterIndex] = uigetfile('*.mat','Select the Mat-file','MultiSelect','on');
% % cc=0;
% % for fcnt=1:35%nsims
% %     cc=cc+1;
% % a=load([PathName,FileName{cc}]);
% % connmat(:,:,:,:,fcnt)=a.parc4_connmat;%(:,:,:,:,fcnt)
% % end
% % end
% % % save('connmat_35sims.mat','connmat')
% %




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