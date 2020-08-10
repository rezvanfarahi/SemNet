clear
clc
close all
%%aparc
%SNR3-chunck1
nsims=36;%length(FileName);
TP_grand=zeros(4,3,5,4,nsims);%the 9: COH:leak3,noise3,leak1,noise1, imCOH:leak3,noise3,ideal,leak1,noise1,
FP_grand=zeros(4,3,5,4,nsims);
TN_grand=zeros(4,3,5,4,nsims);
FN_grand=zeros(4,3,5,4,nsims);
seedcnt=0;
mcc=0;
mc_gmat=zeros(36,5,12);
for nseeds=[3,5,10,15]
    nseeds
    seedcnt=seedcnt+1;
    concnt=0;
    for nconns=[0.25,0.5,1]
        nconns
        mcc=mcc+1;
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
                    load ('perm_pval_parc_nullb_0.mat');
%                     load ('ideal_pval_parc0.mat');
                case 2
                    nnodes=74;
                    inpath=['/imaging/rf02/parcellation/extended_simulations/funcseeds_avgsims70_offset0_seeds',num2str(nseeds),'_conns',num2str(nconns),'_bothsnrs_extglasser/aparc_mod/'];
                    cd(inpath)
                    load ('perm_pval_parc_nullb_1.mat');
%                     load ('ideal_pval_parc1.mat');
                case 3
                    nnodes=148;
                    inpath=['/imaging/rf02/parcellation/extended_simulations/funcseeds_avgsims70_offset0_seeds',num2str(nseeds),'_conns',num2str(nconns),'_bothsnrs_extglasser/aparc2009/'];
                    cd(inpath)
                    load ('perm_pval_parc_nullb_2.mat');
%                     load ('ideal_pval_parc2.mat');
                case 4
                    nnodes=74;
                    inpath=['/imaging/rf02/parcellation/extended_simulations/funcseeds_avgsims70_offset0_seeds',num2str(nseeds),'_conns',num2str(nconns),'_bothsnrs_extglasser/aparc2009_mod/'];
                    cd(inpath)
                    load ('perm_pval_parc_nullb_3.mat');
%                     load ('ideal_pval_parc3.mat');
                case 5
                    nnodes=70;
                    inpath=['/imaging/rf02/parcellation/extended_simulations/funcseeds_avgsims70_offset0_seeds',num2str(nseeds),'_conns',num2str(nconns),'_bothsnrs_extglasser/RG/'];
                    cd(inpath)
                    load ('perm_pval_parc_nullb_4.mat');
%                     load ('ideal_pval_parc4.mat');
            end
                        
            load(['/imaging/rf02/parcellation/extended_simulations/funcseeds_avgsims70_offset0_seeds',num2str(nseeds),'_conns',num2str(nconns),'_bothsnrs_extglasser/allmissedcons.mat'])
%             missedcons=squeeze(mean(missedcons,3));
            missedcons=missedcons(:,which_parc);
            
            
            %% null
%             switch which_parc
%                 case 1
%                     inpath='/imaging/rf02/parcellation/extended_simulations/null_sims40_offset0_bothsnrs/aparc/';
%                     cd(inpath)
%                     aa=dir ('*parc_connmat*.mat');
%                 case 2
%                     inpath='/imaging/rf02/parcellation/extended_simulations/null_sims40_offset0_bothsnrs/aparc_mod/';
%                     cd(inpath)
%                     aa=dir ('*parc_connmat*.mat');
%                 case 3
%                     inpath='/imaging/rf02/parcellation/extended_simulations/null_sims40_offset0_bothsnrs/aparc2009/';
%                     cd(inpath)
%                     aa=dir ('*parc_connmat*.mat');
%                 case 4
%                     inpath='/imaging/rf02/parcellation/extended_simulations/null_sims40_offset0_bothsnrs/aparc2009_mod/';
%                     cd(inpath)
%                     aa=dir ('*parc_connmat*.mat');
%                 case 5
%                     inpath='/imaging/rf02/parcellation/extended_simulations/null_sims40_offset0_bothsnrs/RG/';
%                     cd(inpath)
%                     aa=dir ('*parc_connmat*.mat');
%             end
%             
%             nullconnmat=zeros(nnodes,nnodes,8,nsubs,nsims);
%             %                 [FileName,PathName,FilterIndex] = uigetfile('*.mat',['Select parc',num2str(which_parc),'null'],'MultiSelect','on');
%             cc=0;
%             for fcnt=1:35%nsims
%                 cc=cc+1;
%                 a=load([inpath,aa(fcnt).name]);%load([PathName,FileName{cc}]);
%                 switch which_parc
%                     case 1
%                         nullconnmat(:,:,:,:,fcnt)=a.parc0_connmat;%(:,:,:,:,fcnt)
%                         clear a
%                     case 2
%                         nullconnmat(:,:,:,:,fcnt)=a.parc1_connmat;
%                         clear a
%                     case 3
%                         nullconnmat(:,:,:,:,fcnt)=a.parc2_connmat;
%                         clear a
%                     case 4
%                         nullconnmat(:,:,:,:,fcnt)=a.parc3_connmat;
%                         clear a
%                     case 5
%                         nullconnmat(:,:,:,:,fcnt)=a.parc4_connmat;
%                         clear a
%                 end
%             end
%             nullconnmat=mean(abs(nullconnmat),5);
%             base_connmat=repmat(nullconnmat(:,:,[2,1,1,4,3,6,5,5,8,7],:),[1,1,1,1,nsims]);
%             clear nullmat
            %% null done
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
            disp('significant conns')
%             pfdrmat=zeros(size(pval));
%             for cc1=1:10
% %                 cc1
%                 for cc2=1:nsims
% %                     cc2
%                     pp=pval(:,:,cc1,cc2)+triu(1000*ones(nnodes,nnodes));
%                     pp=pp(pp<1000);
%                     pfdr=mafdr(pp);
%                     ii=0;
%                     for cc3=1:nnodes
%                         for cc4=1:nnodes
%                             if cc3<cc4
%                                 ii=ii+1;
%                                 pfdrmat(cc4,cc3,cc1,cc2)=pfdr(ii);
%                             end
%                         end
%                     end
%                 end
%             end
            ntests=nnodes*(nnodes-1)/2;
            pfdrmat=pval;%*ntests;
            % pthresh=pval/ntests;
            
            
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
%             bcm=squeeze(mean(basecorr_connmat(:,:,3,:,:),4));
            for cc1=1:9
                %                     cc1
                for cc2=1:nsims
                    %                         cc2
                    thismat=squeeze(pfdrmat(:,:,cc1,cc2));
%                     qq=bcm(:,:,cc2)+triu(1000*ones(nnodes,nnodes));
                    thismat(thismat>0.05)=-1;
                    thismat(thismat>-1)=0;
                    if cc1==2 || cc1==5 || cc1==7 || cc1==9 %ideal cases
                        thismat=thismat+1;
                    else
                        thismat=(thismat+2)*2;
                    end
                    thismat=tril(thismat,-1);
                    pfdrmat_final(:,:,cc1,cc2)=thismat;
                end
            end
%             pfdrmat_final(:,:,3,:)=pvalid;
%             eval_mat=pfdrmat_final(:,:,[1,2,4,5,6,7,8,9],:)- pfdrmat_final(:,:,[3,3,3,3,3,3,3,3],:);
%             TP=zeros(8,nsims);
%             TN=zeros(8,nsims);
%             FP=zeros(8,nsims);
%             FN=zeros(8,nsims);
%             for cc1=1:8
%                 for cc2=1:nsims
%                     thiseval=eval_mat(:,:,cc1,cc2);
%                     refvect=pfdrmat_final(:,:,3,cc2)+triu(1000*ones(nnodes,nnodes));
%                     msdcn=missedcons(cc2);
%                     refpos=max([sum(sum(refvect==1)),1]);%+msdcn;
%                     refneg=max([sum(sum(refvect==0)),1]);%-msdcn;
%                     TP(cc1,cc2)=sum(sum(thiseval==3))/refpos;
%                     FP(cc1,cc2)=sum(sum(thiseval==4));%/refneg;%/refpos;
%                     TN(cc1,cc2)=sum(sum(thiseval==4))/refneg;%(sum(sum(thiseval==2))-msdcn)/refneg;
%                     FN(cc1,cc2)=refneg;%(sum(sum(thiseval==1))+msdcn)/refpos;
%                 end
%             end
            eval_mat=pfdrmat_final(:,:,[1,4,6,8],:)- pfdrmat_final(:,:,[2,5,2,5],:);
            TP=zeros(4,nsims);
            TN=zeros(4,nsims);
            FP=zeros(4,nsims);
            FN=zeros(4,nsims);
            for cc1=1:4
                refs=[2,5,2,5];
                for cc2=1:nsims
                    thiseval=eval_mat(:,:,cc1,cc2);
                    refvect=pfdrmat_final(:,:,refs(cc1),cc2)+triu(1000*ones(nnodes,nnodes));
                    msdcn=missedcons(cc2);
                    refpos=max([sum(sum(refvect==1))+msdcn,1]);
                    mc_gmat(cc2,which_parc,mcc)=msdcn/refpos;
                    refneg=max([sum(sum(refvect==0)),1]);%-msdcn;
                    TP(cc1,cc2)=sum(sum(thiseval==3))/refpos;
                    FP(cc1,cc2)=sum(sum(thiseval==4));%/refneg;%/refpos;
                    TN(cc1,cc2)=sum(sum(thiseval==4))/refneg;%(sum(sum(thiseval==2))-msdcn)/refneg;
                    FN(cc1,cc2)=refneg;%(sum(sum(thiseval==1))+msdcn)/refpos;
                end
            end
            
            TP_grand(seedcnt,concnt,which_parc,:,:)=TP;%mean(TP,2); %the 9: COH:leak3,noise3,leak1,noise1, imCOH:leak3,noise3,ideal,leak1,noise1,
            FP_grand(seedcnt,concnt,which_parc,:,:)=FP;%mean(FP,2);
            TN_grand(seedcnt,concnt,which_parc,:,:)=TN;%mean(TN,2);
            FN_grand(seedcnt,concnt,which_parc,:,:)=FN;%mean(FN,2);
            
            
        end
%         output_path='/imaging/rf02/parcellation/extended_simulations/TP_grand.mat';
%         save(output_path,'TP_grand')
%         output_path='/imaging/rf02/parcellation/extended_simulations/TN_grand.mat';
%         save(output_path,'TN_grand')
%         output_path='/imaging/rf02/parcellation/extended_simulations/FP_grand.mat';
%         save(output_path,'FP_grand')
%         output_path='/imaging/rf02/parcellation/extended_simulations/FN_grand.mat';
%         save(output_path,'FN_grand')
    end
end
% figure,
% for cc1=1:size(TP_grand,3)
%         tp=squeeze(TP_grand(:,:,cc1,3));%leak3:1,5/leak1:3,8
%         figure, %subplot(1,size(TP_grand,3),cc1), 
%         imagesc(tp,[0,1]), colormap(hot),colorbar
% end

close all,
% % % wp=3;%leak3:1,5/leak1:3,8 %FP4x3x5x9
% % % for cc1=1:size(TP_grand,2)
% % %         tp1=squeeze(TP_grand(:,cc1,[2,4,5],wp))./squeeze(TP_grand(:,cc1,[1,1,1],wp));%leak3:1,5/leak1:3,8 %FP4x3x5x9
% % %         tp2=squeeze(TP_grand(:,cc1,[2,4,5],wp))./squeeze(TP_grand(:,cc1,[3,3,3],wp));
% % %         figure,%subplot(1,size(TP_grand,2),cc1), 
% % %         imagesc(tp1, [0 max(tp1(:))]), colormap(hot),colorbar
% % %         figure,%subplot(1,size(TP_grand,2),cc1), 
% % %         imagesc(tp2,[0 max(tp1(:))]), colormap(hot),colorbar
% % % end
% % % 
% % % % figure,
% % % % for cc1=1:size(FP_grand,2)
% % % %         fp=squeeze(FP_grand(:,cc1,[1,3,2,4,5],1));%leak3:1,5/leak1:3,8 %FP4x3x5x9
% % % %         figure, %subplot(1,size(FP_grand,2),cc1), 
% % % %         imagesc(fp,[1,100]), colormap(summer),colorbar
% % % % end
% % % 
% % % wp=1;%leak3:1,5/leak1:3,8 %FP4x3x5x9
% % % for cc1=1:size(TP_grand,2)
% % %         tp1=squeeze(FP_grand(:,cc1,[1,1,1],wp))./squeeze(FP_grand(:,cc1,[2,4,5],wp));%leak3:1,5/leak1:3,8 %FP4x3x5x9
% % %         tp2=squeeze(FP_grand(:,cc1,[3,3,3],wp))./squeeze(FP_grand(:,cc1,[2,4,5],wp));
% % %         h(1)=figure('units','normalized','position',[0 0 1 1]);%subplot(1,size(TP_grand,2),cc1), 
% % %         imagesc(tp1,[0 max(tp1(:))]), colormap(summer),h(2)=colorbar;set(h(2),'fontsize',28)
% % %         set(gca,'XTickLabelMode','manual','XTickLabel',[],'XTick',[])
% % %         set(gca,'YTickLabelMode','manual','YTickLabel',[],'YTick',[])
% % % %         set(h(1),'Position',[])
% % %         output_path1='/imaging/rf02/parcellation/extended_simulations/test';
% % %         saveas(gcf,output_path1,'jpg')
% % %         figure,%subplot(1,size(TP_grand,2),cc1), 
% % %         imagesc(tp2,[0 max(tp2(:))]), colormap(summer),colorbar
% % % end
% figure, 
% figure1 = figure('XVisual',...
%     '0x21 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)');
figure1=figure;
% Create axes
axes1 = axes('Parent',figure1,...
    'XTickLabel',{'0','25%','50%','100%','25%','50%','100%','25%','50%','100%','25%','50%','100%','13'},...
    'XTick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13]);
hold(axes1,'all');
cols=[[0 0 1];[1 0 0];[0 0 0];[0.6 0.6 0];[1 0 1]];
marks={'o';'square';'.';'pentagram';'diamond'};
nms={'DKA','mod DKA','DA','mod DA','RG'};
% set(plot1(1),'Marker','o','Color',[0 0 1],'DisplayName','DKA');
% set(plot1(2),'Marker','square','Color',[1 0 0],'DisplayName','mod DKA');
% set(plot1(3),'Marker','^','DisplayName','DA','Color',[0 0 0]);
% set(plot1(4),'Marker','pentagram','Color',[0.6 0.6 0],...
%     'DisplayName','mod DA');
% set(plot1(5),'Marker','diamond','Color',[1 0 1],'DisplayName','RG');
for whichparc=1:5
aa=squeeze(TP_grand(:,:,whichparc,1,:));%./squeeze(TP_grand(:,:,3,3));
aa=permute(aa,[2,1,3]);
bb=reshape(aa,12,nsims);
% bb=bb([1,2,3,4,5,6,7,8,9,10,11,12],:);
bbmean=mean(bb,2);
bbstd=std(bb,0,2);
bbf=[bbmean'-bbstd';bbmean';bbmean'+bbstd'];

aa1=squeeze(TP_grand(:,:,whichparc,1,:));%./squeeze(TP_grand(:,:,3,3));
aa1=permute(aa1,[2,1,3]);
bb1=reshape(aa1,12,nsims);
% bb=bb([1,2,3,4,5,6,7,8,9,10,11,12],:);
bbmean1=mean(bb1,2);
bbstd=std(bb,0,2);
bbf=[bbmean'-bbstd';bbmean';bbmean'+bbstd'];

dpmean=norminv(bbmean)-norminv(bbmean1);
% bb=bb([1,2,4,3,5,6,7,8,10,9,11,12]);
x=1:12;
% hold on, plotshaded(x,bbf,cols(whichparc))%imagesc(bb',[-1 max(bb)]),colormap(hot)
hold on, plot(bbmean,'Color',cols(whichparc,:),'DisplayName',nms{whichparc},'Marker',marks{whichparc})%imagesc(bb',[-1 max(bb)]),colormap(hot)

% hold on, plot(bb,['--',cols(whichparc)])%imagesc(bb',[-1 max(bb)]),colormap(hot)
end
legend(axes1,'show');
cd('/imaging/rf02/parcellation/extended_simulations/results/tpfp')
% p = kruskalwallis([bb(5,:)',bb1(5,:)'])
% pfdrmat(pfdrmat>=0.05)=-13
% pfdrmat(pfdrmat>-1)=0;
% pfdrmat0=pfdrmat+1;
%
% pvalf=zeros(4,5,5);
% tncoh3=squeeze(TN_grand(:,:,:,1,:));
% tncoh3=squeeze(mean(tncoh3,2));
% pvalf=ones(5,5,4);
% for ii=1:5
% for jj=1:5
% for kk=1:4
%     if ii~=jj
% pvalf(ii,jj,kk)=kruskalwallis([squeeze(tncoh3(kk,ii,:)),squeeze(tncoh3(kk,jj,:))],[],'off');
%     end
% close all
% end
% end
% end

tncoh3=squeeze(TP_grand(:,:,:,3,:));
tncoh3=squeeze(mean(tncoh3,2));
tncoh3=squeeze(mean(tncoh3,1));
tncoh3=tncoh3([1,3,2,4,5],:);
mnstd=[mean(tncoh3,2),std(tncoh3,[],2)];
pvalf=ones(5,5);
for ii=1:5
for jj=1:5

    if ii~=jj
pvalf(ii,jj)=kruskalwallis([squeeze(tncoh3(ii,:))',squeeze(tncoh3(jj,:))'],[],'off');
    end
close all

end
end

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