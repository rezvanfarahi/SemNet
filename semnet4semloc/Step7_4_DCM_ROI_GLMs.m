
exp1  = load('/imaging/rf02/Semnet/semnet4semloc/dcm/SemLoc_SD_80wins.mat');
exp2a = load('/imaging/rf02/Semnet/semnet4semloc/dcm/SemNet_SD_80wins.mat');
exp2b = load('/imaging/rf02/Semnet/semnet4semloc/dcm/SemNet_LD_80wins.mat');

Nr = size(exp1.Datafm_wins,1)
Nt = size(exp1.Datafm_wins,2)
xaxis=[50:5:50+79*5];
Ns(1) = size(exp1.Datafm_wins,4);
Ns(2) = size(exp2a.Datafm_wins,4)
Nsubj=Ns(1)+2*Ns(2);

X = [ones(Ns(1),1) zeros(Ns(1),2)];
X = [X; zeros(2*Ns(2),1) kron(eye(2),ones(Ns(2),1))]
X = [X [zeros(Ns(1),Ns(2)); eye(Ns(2)); eye(Ns(2))]];
figure,imagesc(X)

nperms=1000;

Tcom = zeros(Nr,Nt); Fcom = zeros(Nr,Nt); Pcom = zeros(Nr,Nt); 
Tdif = zeros(Nr,Nt); Fdif = zeros(Nr,Nt); Pdif = zeros(Nr,Nt);
Tcom_perm = zeros(Nr,Nt,nperms); Fcom_perm = zeros(Nr,Nt,nperms); Pcom_perm = zeros(Nr,Nt,nperms); 
Tdif_perm = zeros(Nr,Nt,nperms); Fdif_perm = zeros(Nr,Nt,nperms); Pdif_perm = zeros(Nr,Nt,nperms);
load('/imaging/rf02/Semnet/semnet4semloc/dcm/Pdif_80wins.mat')
load('/imaging/rf02/Semnet/semnet4semloc/dcm/Pdif_perm_80wins.mat')
load('/imaging/rf02/Semnet/semnet4semloc/dcm/Fdif_80wins.mat')
load('/imaging/rf02/Semnet/semnet4semloc/dcm/Fdif_perm_80wins.mat')

load('/imaging/rf02/Semnet/semnet4semloc/dcm/Pcom_80wins.mat')
load('/imaging/rf02/Semnet/semnet4semloc/dcm/Pcom_perm_80wins.mat')
load('/imaging/rf02/Semnet/semnet4semloc/dcm/Fcom_80wins.mat')
load('/imaging/rf02/Semnet/semnet4semloc/dcm/Fcom_perm_80wins.mat')
load('/imaging/rf02/Semnet/semnet4semloc/dcm/Tcom_80wins.mat')
load('/imaging/rf02/Semnet/semnet4semloc/dcm/Tcom_perm_80wins.mat')


for r = 1:Nr
    %since we want to cluster over time points,we use the same randomisation order for all of them
    perm_conds=round(rand(nperms,Nsubj));
    perm_conds(perm_conds<1)=-1;
%     perm_conds=zeros(nperms,Nsubj);
%     for pcnt=1:nperms
%         perm_cond=rand(1,Nsubj)-0.5;perm_cond(perm_cond>0)=1;perm_cond(perm_cond<=0)=-1;
%         perm_conds(pcnt,:)=perm_cond; 
%     end
    
    for t = 1:Nt
        y =     squeeze(exp1.Datafm_wins(r,t,:,:))'  * [1 -1]'; % contrast of conditions
        y = [y; squeeze(exp2a.Datafm_wins(r,t,:,:))' * [1 -1]']; % contrast of conditions
        y = [y; squeeze(exp2b.Datafm_wins(r,t,:,:))' * [1 -1]']; % contrast of conditions
        %[Tval,Fcom(r,t),Pcom(r,t),df,R2,cR2,B,res,aR2,iR2,Bcov] = glm(y,X,[1 1 1 zeros(1,Ns(2))]');
        [Tcom(r,t),Fcom(r,t),Pcom(r,t),df,R2,cR2,B,res,aR2,iR2,Bcov] = glm(y,X,[1 1 1 ones(1,Ns(2))*2/Ns(2)]'); % http://www.sbirc.ed.ac.uk/cyril/download/Contrast_Weighting_Glascher_Gitelman_2008.pdf
        [Tdif,Fdif(r,t),Pdif(r,t),df,R2,cR2,B,res,aR2,iR2,Bcov] = glm(y,X,[detrend(eye(3),0) zeros(3,Ns(2))]');        
        for pcnt=1:nperms
            y_perm =y.* perm_conds(pcnt,:)';            
            [Tcom_perm(r,t,pcnt),Fcom_perm(r,t,pcnt),Pcom_perm(r,t,pcnt),df,R2,cR2,B,res,aR2,iR2,Bcov] = glm(y_perm,X,[1 1 1 ones(1,Ns(2))*2/Ns(2)]',-1); % http://www.sbirc.ed.ac.uk/cyril/download/Contrast_Weighting_Glascher_Gitelman_2008.pdf
            [Tdif_perm,Fdif_perm(r,t,pcnt),Pdif_perm(r,t,pcnt),df,R2,cR2,B,res,aR2,iR2,Bcov] = glm(y_perm,X,[detrend(eye(3),0) zeros(3,Ns(2))]',-1);
        end
    end
end

alpha = 0.05 % uncorrected
%alpha = alpha/(Nr*Nt) % extreme Bonferonni, ignoring that time windows correlated
%alpha = alpha/Nr % compromise!

figure,imagesc(Pcom),colorbar,xlabel('Twin'),ylabel('ROI'),set(gca,'XTick',[1:10:Nt],'XTickLabel',[50:50:445],'YTick',[1:5]),colormap('gray'),caxis([0 alpha]),title('Common condition effect across exps')
figure,imagesc(Pdif),colorbar,xlabel('Twin'),ylabel('ROI'),set(gca,'XTick',[1:10:Nt],'XTickLabel',[50:50:445],'YTick',[1:5]),colormap('gray'),caxis([0 alpha]),title('Interaction between exp and condition')

[cpvals_com,cstart_com,cend_com,croi_com]=cluster_perm(Pcom,Pcom_perm,Tcom,Tcom_perm,xaxis);
[cpvals_dif,cstart_dif,cend_dif,croi_dif]=cluster_perm(Pdif,Pdif_perm,Fdif,Fdif_perm,xaxis);

% % % %% cluster-based permutation
% % % Pdif_bin=Pdif;
% % % Pdif_bin(Pdif_bin>=0.05)=1;Pdif_bin(Pdif_bin<1)=0;Pdif_bin=1-Pdif_bin;
% % % Pdif_perm_bin=Pdif_perm;
% % % Pdif_perm_bin(Pdif_perm_bin>=0.05)=1;Pdif_perm_bin(Pdif_perm_bin<1)=0;Pdif_perm_bin=1-Pdif_perm_bin;
% % % for r = 1:Nr
% % %     pdif_row=Pdif_bin(r,:);
% % %     tdif_row=Fdif(r,:);
% % %      
% % %     %find clusters in real data
% % %     clu = pdif_row ~= 0;
% % %     clu_b=findstr([0 clu], [0 1]);
% % %     clu_e=findstr([clu 0], [1 0]);
% % %     tcluster=[];
% % %     for cc=1:length(clu_b)
% % %         tcluster(cc)=sum(tdif_row(clu_b(cc):clu_e(cc)));
% % %     end
% % %     %find largest cluster in surrogate data
% % %     tcluster_perm=[];
% % %     for pcnt=1:nperms
% % %         pdif_row=squeeze(Pdif_perm_bin(r,:,pcnt));
% % %         tdif_row=squeeze(Fdif_perm(r,:,pcnt));
% % %         
% % %         clu = pdif_row ~= 0;
% % %         if sum(clu)>0
% % %             clu_bp=findstr([0 clu], [0 1]);
% % %             clu_ep=findstr([clu 0], [1 0]);
% % %             tclu_perm=[];
% % %             for cc=1:length(clu_bp)
% % %                 tclu_perm(cc)=sum(tdif_row(clu_bp(cc):clu_ep(cc)));
% % %             end
% % %             tcluster_perm(pcnt)=max(abs(tclu_perm));
% % %         else
% % %             tcluster_perm(pcnt)=0;
% % %         end
% % %         
% % %     end
% % %     cluster_pvalue=[];
% % %     for cnt=1:length(tcluster)
% % %         cluster_pvalue=sum(abs(tcluster(cnt))<tcluster_perm)/nperms;
% % %         if cluster_pvalue<0.05
% % %             fprintf('cluster-pvalue=%f, start=%d,  end=%d, ROI=%d \n',cluster_pvalue,xaxis(clu_b(cnt)),xaxis(clu_e(cnt)),r);       
% % %         end
% % %     end
% % %     
% % %     
% % %     
% % % end
