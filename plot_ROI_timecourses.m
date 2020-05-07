clear
clc
close all
% bands={'alpha','gamma'};
% for bcnt=1:length(bands)
%     band=bands{bcnt};
% load(['/imaging/rf02/Semnet/ROI_analysis/',band,'_newfilt_wbw_ROItcs_normal_meanflip_ROI_ATAGIFGMTGVISHRHN_evVisHrHn.mat'])
% tstep=1;
% aa=1:tstep:size(tcmat_normf,3)-tstep;
% bb=1+tstep:tstep:size(tcmat_normf,3);
% pmm=zeros(20,7,length(aa),3);
% for ii=1:length(aa)
% pmm(:,:,ii,:)=squeeze(mean(tcmat_normf(:,:,aa(ii):bb(ii),:),3));
% end
% pm=squeeze(mean(pmm,1));
% tline=-200:tstep:450-tstep;
% tvect=100/tstep:749/tstep;
% figure, plot(tline,squeeze(pm(5,tvect,1)),'k')
% hold on, plot(tline,squeeze(pm(5,tvect,2)),'b')
% hold on, plot(tline,squeeze(pm(5,tvect,3)),'r')
% title([band,' visual cortex'])
% figure, plot(tline,squeeze(pm(7,tvect,1)),'k')
% hold on, plot(tline,squeeze(pm(7,tvect,2)),'b')
% hold on, plot(tline,squeeze(pm(7,tvect,3)),'r')
% title([band,' motor cortex'])
% figure, plot(tline,squeeze(pm(6,tvect,1)),'k')
% hold on, plot(tline,squeeze(pm(6,tvect,2)),'b')
% hold on, plot(tline,squeeze(pm(6,tvect,3)),'r')
% title([band,' auditory cortex'])
% end
% 
% 
% close all
%%broad band
load(['/imaging/rf02/Semnet/ROI_analysis/wbw_ROItcs_normal_meanflip_to30_ROI_ATAGIFGMTGVISHRHN_evVisHrHn3.mat'])
nsubs=19;
nrois=7;
tstep=1;
aa=1:tstep:size(tcmat_normf,3)-tstep;
bb=1+tstep:tstep:size(tcmat_normf,3);
pmm=zeros(nsubs,nrois,length(aa),3);
for ii=1:length(aa)
pmm(:,:,ii,:)=squeeze(mean(tcmat_normf(:,:,aa(ii):bb(ii),:),3));
end
pm=squeeze(mean(pmm,1));
tline=-200:tstep:450-tstep;
tvect=100/tstep:749/tstep;
figure, plot(tline,squeeze(pm(5,tvect,1)),'r')
hold on, plot(tline,squeeze(pm(5,tvect,3)),'k')
hold on, plot(tline,squeeze(pm(7,tvect,1)),'r')
hold on, plot(tline,squeeze(pm(7,tvect,3)),'k')
title(['Visual-Hand Condition x ROI Interaction Effect'])
const=nan(size(tvect));
const(400:475)=-3.5*(10^-12);
hold on, plot(tline,const,'b')
const1=nan(size(tvect));
const1(350:375)=-3.5*(10^-12);
hold on, plot(tline,const1,'c.')

figure, plot(tline,squeeze(pm(5,tvect,1)),'r')
hold on, plot(tline,squeeze(pm(5,tvect,2)),'b')
hold on, plot(tline,squeeze(pm(6,tvect,1)),'r')
hold on, plot(tline,squeeze(pm(6,tvect,2)),'b')
title(['Visual-Auditory Condition x ROI Interaction Effect'])
const1=nan(size(tvect));
const1(525:550)=-3.5*(10^-12);
hold on, plot(tline,const1,'b')

figure, plot(tline,squeeze(pm(6,tvect,2)),'b')
hold on, plot(tline,squeeze(pm(6,tvect,3)),'k')
hold on, plot(tline,squeeze(pm(7,tvect,2)),'b')
hold on, plot(tline,squeeze(pm(7,tvect,3)),'k')
title(['Hand-Auditory Condition x ROI Interaction Effect'])

% figure, plot(tline,squeeze(pm(5,tvect,1)),'r')
% hold on, plot(tline,squeeze(pm(5,tvect,2)),'b')
% hold on, plot(tline,squeeze(pm(5,tvect,3)),'k')
% title('visual cortex')
% figure, plot(tline,squeeze(pm(7,tvect,1)),'r')
% hold on, plot(tline,squeeze(pm(7,tvect,2)),'b')
% hold on, plot(tline,squeeze(pm(7,tvect,3)),'k')
% title(['motor cortex'])
% figure, plot(tline,squeeze(pm(6,tvect,1)),'r')
% hold on, plot(tline,squeeze(pm(6,tvect,2)),'b')
% hold on, plot(tline,squeeze(pm(6,tvect,3)),'k')
% title(['auditory cortex'])


