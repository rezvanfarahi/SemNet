
clear
clc
close all
addpath(genpath("/Users/rezvanh/Documents/PhD/boundedline_pkg"))
addpath(genpath("/Users/rezvanh/Documents/PhD/shadedErrorBar"))
%ref for shaded area: https://uk.mathworks.com/matlabcentral/answers/180829-shade-area-between-graphs

%refs for standard error bar: 
%1 within subject error bar; easy ref: http://www.cogsci.nl/blog/tutorials/156-an-easy-way-to-create-graphs-with-within-subject-error-bars
%main ref: http://www.tqmp.org/RegularArticles/vol01-1/p042/p042.pdf

%2 cross-subject error bar; var1 corresponding to semloc, 
%then I treated task2,3 as one (same semnet subjs) and found
%one var23 for them, then computed pooled variance: https://en.wikipedia.org/wiki/Pooled_variance
%also useful for pooled var- answer2: https://stats.stackexchange.com/questions/231027/combining-samples-based-off-mean-and-standard-error
%combined standard error is based on pooled variance

inpath='/Users/rezvanh/Documents/PhD/manuscripts/semloc_july2020/semnet4semloc/';
filenames={'SemLoc_SD_500wins.mat','SemNet_SD_500wins.mat','SemNet_LD_500wins.mat'};
nroi=5;
ntime=500;
ncond=2;
nsubj=53;
lband=34;
hband=68;
X=zeros(nroi,ntime,ncond,nsubj);
for ii=1:3
    inname=[inpath,filenames{ii}];
    load(inname)
    if ii==1
        X(:,:,:,1:17)=Datafm_wins;
    end
    if ii==2
        X(:,:,:,18:35)=Datafm_wins;
    end
    if ii==3
        X(:,:,:,36:end)=Datafm_wins;
    end
    
end
X=abs(X);
Xm=mean(X(:,:,:,:),4);
X25=zeros(size(Xm));
X75=zeros(size(Xm));
Xerr=zeros(size(Xm));
nsubj1=17;
nsubj2=(nsubj-17)/2;
nsubj3=(nsubj-17)/2;

for ii=1:nroi
    for jj=1:ntime
        subavg1=squeeze(mean(X(ii,jj,:,1:17),3));gavg1=mean(subavg1);
        subavg2=squeeze(mean(X(ii,jj,:,18:35),3));gavg2=mean(subavg2);
        subavg3=squeeze(mean(X(ii,jj,:,36:end),3));gavg3=mean(subavg3);
        subavg=squeeze(mean(X(ii,jj,:,:),3));gavg=mean(subavg);
        for kk=1:ncond
            X25(ii,jj,kk)=prctile(squeeze(X(ii,jj,kk,:)),lband);
            X75(ii,jj,kk)=prctile(squeeze(X(ii,jj,kk,:)),hband);
            nval1=squeeze(X(ii,jj,kk,1:17))-subavg1+gavg1;
            nval2=squeeze(X(ii,jj,kk,18:35))-subavg2+gavg2;
            nval3=squeeze(X(ii,jj,kk,36:end))-subavg3+gavg3;
            nval23=[nval2;nval3];
            v1=(nsubj1-1)*var(nval1);
            v2=(nsubj2-1)*var(nval2);
            v3=(nsubj3-1)*var(nval3);
            v23=(nsubj2+nsubj3-1)*var(nval23);
            varpool=(v1+v23)/(nsubj-2);
            Xerr(ii,jj,kk)=sqrt(varpool/nsubj);%1.96*
%             if kk<18
%                 Xerr(ii,jj,kk)=std(squeeze(X(ii,jj,kk,:))-subavg+gavg1)/sqrt(nsubj);
%             elseif kk<36
%                 Xerr(ii,jj,kk)=std(squeeze(X(ii,jj,kk,:))-subavg+gavg2)/sqrt(nsubj);
%             else
%                 Xerr(ii,jj,kk)=std(squeeze(X(ii,jj,kk,:))-subavg+gavg3)/sqrt(nsubj);
%             end
        end
    end
end
%% plot ERP- no error bar
for ii=1:5
    figure();
    conmat=squeeze(Xm(ii,:,1));
    absmat=squeeze(Xm(ii,:,2));
    plot(-50:1:449,conmat,'r','LineWidth',5)
    hold on   
    plot(-50:1:449,absmat,'b','LineWidth',5)
    set(gca,'FontSize',20)
    legend('Concrete','Abstract')
end

%% plot with error bars
tmin=50;
for ii=1:5
    figure();
    conmean=squeeze(Xm(ii,tmin:end,1));
    %constd=squeeze(Xstd(ii,:,1));
    con25=conmean-squeeze(Xerr(ii,tmin:end,1));
    con75=conmean+squeeze(Xerr(ii,tmin:end,1));
    
    absmean=squeeze(Xm(ii,tmin:end,2));
    %absstd=squeeze(Xstd(ii,:,2));
    abs25=absmean-squeeze(Xerr(ii,tmin:end,2));
    abs75=absmean+squeeze(Xerr(ii,tmin:end,2));
    x=-50+tmin-1:1:449;
    plot(x, con25, '--r', 'LineWidth', 1);
    hold on;
    plot(x, con75, '--r', 'LineWidth', 1);
    x2 = [x, fliplr(x)];
    inBetween = [con25, fliplr(con75)];
    fill(x2, inBetween, 'r','facealpha',0.2);
    plot(x,conmean,'r','LineWidth', 2)
    
    plot(x, abs25, '--b', 'LineWidth', 1);
    hold on;
    plot(x, abs75, '--b', 'LineWidth', 1);
    x2 = [x, fliplr(x)];
    inBetween = [abs25, fliplr(abs75)];
    fill(x2, inBetween, 'b','facealpha',0.2);
    plot(x,absmean,'b','LineWidth', 2)
    set(gca,'FontSize',20)
    %legend('Concrete','Abstract')
 
end
%% plot with error bars based on percentile- useless I think
% for ii=1:5
%     figure();
%     conmean=squeeze(Xm(ii,:,1));
%     %constd=squeeze(Xstd(ii,:,1));
%     con25=squeeze(X25(ii,:,1));
%     con75=squeeze(X75(ii,:,1));
%     conperc=[con25;con75]';
%     absmean=squeeze(Xm(ii,:,2));
%     %absstd=squeeze(Xstd(ii,:,2));
%     abs25=squeeze(X25(ii,:,2));
%     abs75=squeeze(X75(ii,:,2));
%     absperc=[abs25;abs75]';
%     x=-50:1:449;
%     plot(x, con25, '--r', 'LineWidth', 1);
%     hold on;
%     plot(x, con75, '--r', 'LineWidth', 1);
%     x2 = [x, fliplr(x)];
%     inBetween = [con25, fliplr(con75)];
%     fill(x2, inBetween, 'r','facealpha',0.2);
%     plot(x,conmean,'r','LineWidth', 2)
%     
%     plot(x, abs25, '--b', 'LineWidth', 1);
%     hold on;
%     plot(x, abs75, '--b', 'LineWidth', 1);
%     x2 = [x, fliplr(x)];
%     inBetween = [abs25, fliplr(abs75)];
%     fill(x2, inBetween, 'b','facealpha',0.2);
%     plot(x,absmean,'b','LineWidth', 2)
%  
% end
cd '/Users/rezvanh/Documents/PhD/manuscripts/semloc_july2020/semnet4semloc/jpg/evoked_ROI'

   
    
    


    %errorbar(x,conmean,constd)
    %shadedErrorBar(x, conmean, conperc');%, 'r', x, absmean, absperc, 'b');
%     plot(x,conmean,'r')
%     inBetween = [con25', con75'];
%     x2 = [x', fliplr(x)'];   
%     h=fill(con25, con75, 'r','edgecolor','none');
%     set(h,'facealpha',.5)
    
%     curve1 = log(x);
%     curve2 = 2*log(x);
%     plot(x, curve1, 'r', 'LineWidth', 2);
%     hold on;
%     plot(x, curve2, 'b', 'LineWidth', 2);
%     
%     inBetween = [curve1, fliplr(curve2)];
%     