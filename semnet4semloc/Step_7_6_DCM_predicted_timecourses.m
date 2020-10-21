
clear
clc
close all
inpath='/Users/rezvanh/Documents/PhD/manuscripts/semloc_july2020/semnet4semloc/';
filenames={'SemLoc_SD_predicted_250ms_ATL.mat','SemNet_SD_predicted_250ms_ATL.mat','SemNet_LD_predicted_250ms_ATL.mat'};
%filenames={'SemLoc_SD_predicted_450ms_AG.mat','SemNet_SD_predicted_450ms_AG.mat','SemNet_LD_predicted_450ms_AG.mat'};

tmax=251;
X=zeros(tmax,5,53,2);
for ii=1:3
    inname=[inpath,filenames{ii}];
    load(inname)
    if ii==1
        X(:,:,1:17,:)=H_CA;
    end
    if ii==2
        X(:,:,18:35,:)=H_CA;
    end
    if ii==3
        X(:,:,36:end,:)=H_CA;
    end
    
end
Xm=squeeze(mean(X,3));
titles={'ATL','IFG','MTG','AG','vWFA'};
for ii=1:5
    figure();
    conmat=squeeze(Xm(:,ii,1));
    absmat=squeeze(Xm(:,ii,2));
    plot(1:1:tmax,conmat,'r','LineWidth',5)
    hold on   
    plot(1:1:tmax,absmat,'b','LineWidth',5)
    set(gca,'FontSize',45)
    %title(titles{ii})
    %legend('Concrete','Abstract')
end
cd '/Users/rezvanh/Documents/PhD/manuscripts/semloc_july2020/semnet4semloc/jpg/evoked_ROI'