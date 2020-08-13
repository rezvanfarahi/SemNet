function [cluster_pvalues,cstart,cend,rois] = cluster_perm(Pval,Pval_perm,Tval,Tval_perm,xaxis);

% a quick implementation of cluster-based permutation over a time series


Pval_bin=Pval;
Pval_bin(Pval_bin>=0.05)=1;Pval_bin(Pval_bin<1)=0;Pval_bin=1-Pval_bin;
Pval_perm_bin=Pval_perm;
Pval_perm_bin(Pval_perm_bin>=0.05)=1;Pval_perm_bin(Pval_perm_bin<1)=0;Pval_perm_bin=1-Pval_perm_bin;
Nr=size(Pval,1);
nperms=size(Pval_perm,3);
cluster_pvalues=[];
cstart=[];
cend=[];
rois=[];
cpcnt=0;
for r = 1:Nr
    Pval_row=Pval_bin(r,:);
    Tval_row=Tval(r,:);
     
    %find clusters in real data
    clu = Pval_row ~= 0;
    clu_b=findstr([0 clu], [0 1]);
    clu_e=findstr([clu 0], [1 0]);
    tcluster=[];
    for cc=1:length(clu_b)
        tcluster(cc)=sum(Tval_row(clu_b(cc):clu_e(cc)));
    end
    %find largest cluster in surrogate data
    tcluster_perm=[];
    for pcnt=1:nperms
        Pval_row=squeeze(Pval_perm_bin(r,:,pcnt));
        Tval_row=squeeze(Tval_perm(r,:,pcnt));
        
        clu = Pval_row ~= 0;
        if sum(clu)>0
            clu_bp=findstr([0 clu], [0 1]);
            clu_ep=findstr([clu 0], [1 0]);
            tclu_perm=[];
            for cc=1:length(clu_bp)
                tclu_perm(cc)=sum(Tval_row(clu_bp(cc):clu_ep(cc)));
            end
            tcluster_perm(pcnt)=max(abs(tclu_perm));
        else
            tcluster_perm(pcnt)=0;
        end
        
    end
    
    for cnt=1:length(tcluster)
        cluster_pvalue=sum(abs(tcluster(cnt))<tcluster_perm)/nperms;
        if cluster_pvalue<0.05
            cpcnt=cpcnt+1;
            cluster_pvalues(cpcnt)=cluster_pvalue;
            cstart(cpcnt)=xaxis(clu_b(cnt));
            cend(cpcnt)=xaxis(clu_e(cnt));
            rois(cpcnt)=r;
            fprintf('cluster-pvalue=%3.3f, start=%d,  end=%d, ROI=%d \n',cluster_pvalue,xaxis(clu_b(cnt)),xaxis(clu_e(cnt)),r);       
        end
    end    
    
end


return

