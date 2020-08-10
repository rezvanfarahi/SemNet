clear 
clc 
%[num,txt,raw] = xlsread('C:\Users\Rezvan\Documents\PhD period\semloc\first year report\first year report\valence\new2013.xlsx');
[num,txt,raw] = xlsread('C:\Users\rf02\Documents\Rezvan\PhDproject\experimentdesign\nextexperiment\new2013.xlsx');%for emoth
%[num,txt,raw] =
%xlsread('C:\Users\rf02\Documents\Rezvan\PhDproject\experimentdesign\nextex
%periment\AoA_51715_words.xlsx'); for AoA

%[num1,txt1,raw1] = xlsread('C:\Users\Rezvan\Downloads\AffectiveNormsforEnglishWords.xls',2);
%[num1,txt1,raw1] = xlsread('C:\Users\Rezvan\Documents\PhD period\experiment design\new experiment\abstract.xlsx');
[num1,txt1,raw1] = xlsread('C:\Users\rf02\Documents\Rezvan\PhDproject\experimentdesign\nextexperiment\new_words_modified_06Nov.xlsx',6);


conabs=raw1(2:end,1);%reshape(txt1,93*2,1);
anewlist=raw(2:end,2);
anewscale=raw(2:end,[3,4,6,7]);%11 for AoA, [3,4,6,7] for emosh
 for ii=1:size(conabs,1)
ii
for jj =1:length(anewscale)
if strcmpi(conabs(ii,1),anewlist(jj))
conabs(ii,2:5)=anewscale(jj,:);%conabs(ii,2) for AoA, conabs(ii,2:5) for emosh
end
end
end

