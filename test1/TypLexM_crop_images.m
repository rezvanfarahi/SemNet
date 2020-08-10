clear
clc
% cd('/imaging/rf02/TypLexMEG/icaanalysis_results/stc')
% [fn,fp,fi]=uigetfile('*.stc','MultiSelect','on')
% for ii=1:length(fn)/2
% fprintf([fn{ii*2}(1:end-7),'\n'])
% end
% cd('\\cbsu\data\Imaging\rf02\TypLexMEG\icaanalysis_results\jpg\permutation\connectivity\');
%cd('/imaging/rf02/TypLexMEG/icaanalysis_results/jpg/permutation/power/final_07oct15')%connectivity/WB_spokes/spectral/')
% cd('/imaging/rf02/TypLexMEG/icaanalysis_results/jpg/permutation/connectivity/WB_hubs/temporal')%connectivity/WB_spokes/spectral/')
% cd('/imaging/rf02/TypLexMEG/icaanalysis_results/jpg/UVttest/connectivity/WB_hubs/st')
cd('/imaging/rf02/TypLexMEG/icaanalysis_results/typlex/jpg/permutation/evoked')
%cd('/imaging/rf02/TypLexMEG/icaanalysis_results/typlex/jpg/UVttest/power/')

[FileName,PathName] = uigetfile('*.jpg','Select jpg file','MultiSelect', 'on');
for ii=1:length(FileName)
    ii
I = imread([PathName,FileName{ii}]);
if size(I,2)>400
    I2 = imcrop(I,[110 66 380 268]); %%for whole brain images
%     I2 = imcrop(I,[40,120,970,970]);%[120 180 1180 1150]);
%     I2=imcrop(I,[30 40 1160 1100]); %% for circmaps old
%I2=imcrop(I,[175 190 880 825]);
     imwrite(I2,[PathName,FileName{ii}])
end
end