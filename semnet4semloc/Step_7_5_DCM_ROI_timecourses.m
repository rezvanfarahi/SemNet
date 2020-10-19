
clear
clc
close all
inpath='/Users/rezvanh/Documents/PhD/manuscripts/semloc_july2020/semnet4semloc/';
filenames={'SemLoc_SD_500wins.mat','SemNet_SD_500wins.mat','SemNet_LD_500wins.mat'};
X=zeros(5,500,2,53);
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
Xm=mean(X,4);
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
