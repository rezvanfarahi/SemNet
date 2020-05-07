clear
clc
% cd('/imaging/rf02/TypLexMEG/icaanalysis_results/stc')
% [fn,fp,fi]=uigetfile('*.stc','MultiSelect','on')
% for ii=1:length(fn)/2
% fprintf([fn{ii*2}(1:end-7),'\n'])
% end
% cd('\\cbsu\data\Imaging\rf02\TypLexMEG\icaanalysis_results\jpg\permutation\connectivity\');

% cd('/imaging/rf02/Semnet/jpg/permutation/evoked/bands/newregsigned/taugamma_newev/Hear/')
cd('/imaging/rf02/Semnet/jpg/uvttest/evoked/bands/newfilt/newregsigned/vizcats_xband/vishear/')
%cd('/imaging/rf02/TypLexMEG/icaanalysis_results/typlex/jpg/UVttest/power/')
PathNameo='/imaging/rf02/Semnet/jpg/uvttest/evoked/bands/newfilt/newregsigned/vizcats_xband/vishear/cropped/';
[FileName,PathName] = uigetfile('*.jpg','Select jpg file','MultiSelect', 'on');
for ii=1:length(FileName)
    ii
I = imread([PathName,FileName{ii}]);
FileNamet=['cropped',FileName{ii}];

if size(I,2)>400
%     I2 = imcrop(I,[110 66 380 268]); %%for whole brain images
if strcmp(FileName{ii}(11:13),'lat')%(1:3)
    I2 = imcrop(I,[271 156 845 562]); %%for mne created whole brain images
    %[305,667,785 ,537]
else
        I2 = imcrop(I,[271 133 845 566]); %%for mne created whole brain images
        %[281,130]
end

%     I2 = imcrop(I,[40,120,970,970]);%[120 180 1180 1150]);
%     I2=imcrop(I,[30 40 1160 1100]); %% for circmaps old
%I2=imcrop(I,[175 190 880 825]);
     imwrite(I2,[PathNameo,FileNamet])
end
end