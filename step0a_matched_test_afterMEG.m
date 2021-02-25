clear
clc
close all
%% Orig
% [numv,txtv,rawv] =xlsread('C:\Users\rf02\Documents\Rezvan\PhDproject\experimentdesign\nextexperiment\rating\finaldata\matcholder\40\visual40.xlsx');
% [numhr,txthr,rawhr] =xlsread('C:\Users\rf02\Documents\Rezvan\PhDproject\experimentdesign\nextexperiment\rating\finaldata\matcholder\40\hear40.xlsx');
% [numhn,txthn,rawhn] =xlsread('C:\Users\rf02\Documents\Rezvan\PhDproject\experimentdesign\nextexperiment\rating\finaldata\matcholder\40\hand40.xlsx');
% [numna,txtna,rawna] =xlsread('C:\Users\rf02\Documents\Rezvan\PhDproject\experimentdesign\nextexperiment\rating\finaldata\matcholder\40\nabs40.xlsx');
% [numea,txtea,rawea] =xlsread('C:\Users\rf02\Documents\Rezvan\PhDproject\experimentdesign\nextexperiment\rating\finaldata\matcholder\40\eabs40.xlsx');
% labels=[repmat({'visual'},40,1);repmat({'hear'},40,1);repmat({'hand'},40,1);repmat({'neutral'},40,1);repmat({'emotional'},40,1)];
% %word	1 concrete	2 visual	3 hand	4 hear	5 valence	6 AoA	7 LEN	8 FREQ	9 Orth	10 UN2_C	11UN3_C
% %boxplot([numv(:,5);numhr(:,5);numhn(:,5);numna(:,5);numea(:,5)],labels)
% for ii=10
% 
% % anova1([numv(:,ii);numea(:,ii)],labels(1:100))
% anova1([numv(:,ii);numhr(:,ii);numhn(:,ii);numna(:,ii);numea(:,ii)],labels)
% end


%% handhear

% 
% [numhr,txthr,rawhr] =xlsread('C:\Users\rf02\Documents\Rezvan\PhDproject\experimentdesign\nextexperiment\rating\finaldata\matcholder\handhear\hear50.xlsx');
% [numhn,txthn,rawhn] =xlsread('C:\Users\rf02\Documents\Rezvan\PhDproject\experimentdesign\nextexperiment\rating\finaldata\matcholder\handhear\hand50.xlsx');
% [numna,txtna,rawna] =xlsread('C:\Users\rf02\Documents\Rezvan\PhDproject\experimentdesign\nextexperiment\rating\finaldata\matcholder\handhear\nabs50.xlsx');
% [numea,txtea,rawea] =xlsread('C:\Users\rf02\Documents\Rezvan\PhDproject\experimentdesign\nextexperiment\rating\finaldata\matcholder\handhear\eabs50.xlsx');
% labels=[repmat({'hear'},50,1);repmat({'hand'},50,1);repmat({'neutral'},50,1);repmat({'emotional'},50,1)];
% %word	1 concrete	2 visual	3 hand	4 hear	5 valence	6 AoA	7 LEN	8 FREQ	9 Orth	10 UN2_C	11UN3_C
% %boxplot([numv(:,5);numhr(:,5);numhn(:,5);numna(:,5);numea(:,5)],labels)
% for ii=5
% 
% % anova1([numv(:,ii);numea(:,ii)],labels(1:100))
% anova1([numhr(:,ii);numhn(:,ii);numna(:,ii);numea(:,ii)],labels)
% end

% %% text read
% fname=['C:\Users\rf02\Documents\Rezvan\PhDproject\experimentdesign\nextexperiment\rating\finaldata\matcholder\50\word_eabs_matched.txt'];
% [word,concrete,visual,hand,hear,valence,AoA,LEN1,FREQ1,OrthN,UN2_F,UN3_F]=textread(fname,'%s %f %f %f %f %f %f %f %f %f %f %f','headerlines',1);
% numea=[concrete,visual,hand,hear,valence,AoA,LEN1,FREQ1,OrthN,UN2_F,UN3_F];
% fname=['C:\Users\rf02\Documents\Rezvan\PhDproject\experimentdesign\nextexperiment\rating\finaldata\matcholder\50\word_nabs_matched.txt'];
% [word,concrete,visual,hand,hear,valence,AoA,LEN1,FREQ1,OrthN,UN2_F,UN3_F]=textread(fname,'%s %f %f %f %f %f %f %f %f %f %f %f','headerlines',1);
% numna=[concrete,visual,hand,hear,valence,AoA,LEN1,FREQ1,OrthN,UN2_F,UN3_F];
% fname=['C:\Users\rf02\Documents\Rezvan\PhDproject\experimentdesign\nextexperiment\rating\finaldata\matcholder\50\word_vis_matched.txt'];
% [word,concrete,visual,hand,hear,valence,AoA,LEN1,FREQ1,OrthN,UN2_F,UN3_F]=textread(fname,'%s %f %f %f %f %f %f %f %f %f %f %f','headerlines',1);
% numv=[concrete,visual,hand,hear,valence,AoA,LEN1,FREQ1,OrthN,UN2_F,UN3_F];
% fname=['C:\Users\rf02\Documents\Rezvan\PhDproject\experimentdesign\nextexperiment\rating\finaldata\matcholder\50\word_hear_matched.txt'];
% [word,concrete,visual,hand,hear,valence,AoA,LEN1,FREQ1,OrthN,UN2_F,UN3_F]=textread(fname,'%s %f %f %f %f %f %f %f %f %f %f %f','headerlines',1);
% numhr=[concrete,visual,hand,hear,valence,AoA,LEN1,FREQ1,OrthN,UN2_F,UN3_F];
% fname=['C:\Users\rf02\Documents\Rezvan\PhDproject\experimentdesign\nextexperiment\rating\finaldata\matcholder\50\word_hand_matched.txt'];
% [word,concrete,visual,hand,hear,valence,AoA,LEN1,FREQ1,OrthN,UN2_F,UN3_F]=textread(fname,'%s %f %f %f %f %f %f %f %f %f %f %f','headerlines',1);
% numhn=[concrete,visual,hand,hear,valence,AoA,LEN1,FREQ1,OrthN,UN2_F,UN3_F];
% 
% labels=[repmat({'visual'},50,1);repmat({'hear'},50,1);repmat({'hand'},50,1);repmat({'neutral'},50,1);repmat({'emotional'},50,1)];
% %word	1 concrete	2 visual	3 hand	4 hear	5 valence	6 AoA	7 LEN	8 FREQ	9 Orth	10 UN2_C	11UN3_C
% %boxplot([numv(:,5);numhr(:,5);numhn(:,5);numna(:,5);numea(:,5)],labels)
% for ii=7
% 
% % anova1([numv(:,ii);numea(:,ii)],labels(1:100))
% anova1([numv(:,ii);numhr(:,ii);numhn(:,ii);numna(:,ii);numea(:,ii)],labels)
% end

%% text read 3cat

% % fname=['C:\Users\rf02\Documents\Rezvan\PhDproject\experimentdesign\nextexperiment\rating\finaldata\matcholder\3cat\word_vis_matched3.txt'];
% fname=['C:\Users\rf02\Documents\Rezvan\PhDproject\experimentdesign\nextexperiment\MEG\MEG_pilot\wordlist\word_vis2_matched3.txt'];
% [word,concrete,visual,hand,hear,valence,AoA,LEN1,FREQ1,OrthN,UN2_F,UN3_F]=textread(fname,'%s %f %f %f %f %f %f %f %f %f %f %f','headerlines',1);
% numv=[concrete,visual,hand,hear,valence,AoA,LEN1,FREQ1,OrthN,UN2_F,UN3_F];
% % fname=['C:\Users\rf02\Documents\Rezvan\PhDproject\experimentdesign\nextexperiment\rating\finaldata\matcholder\3cat\word_hear_matched3.txt'];
% fname=['C:\Users\rf02\Documents\Rezvan\PhDproject\experimentdesign\nextexperiment\MEG\MEG_pilot\wordlist\word_hear2_matched3.txt'];
% [word,concrete,visual,hand,hear,valence,AoA,LEN1,FREQ1,OrthN,UN2_F,UN3_F]=textread(fname,'%s %f %f %f %f %f %f %f %f %f %f %f','headerlines',1);
% numhr=[concrete,visual,hand,hear,valence,AoA,LEN1,FREQ1,OrthN,UN2_F,UN3_F];
% % fname=['C:\Users\rf02\Documents\Rezvan\PhDproject\experimentdesign\nextexperiment\rating\finaldata\matcholder\3cat\word_hand_matched3.txt'];
% fname=['C:\Users\rf02\Documents\Rezvan\PhDproject\experimentdesign\nextexperiment\MEG\MEG_pilot\wordlist\word_hand2_matched3.txt'];
% [word,concrete,visual,hand,hear,valence,AoA,LEN1,FREQ1,OrthN,UN2_F,UN3_F]=textread(fname,'%s %f %f %f %f %f %f %f %f %f %f %f','headerlines',1);
% numhn=[concrete,visual,hand,hear,valence,AoA,LEN1,FREQ1,OrthN,UN2_F,UN3_F];
% fname=['C:\Users\rf02\Documents\Rezvan\PhDproject\experimentdesign\nextexperiment\rating\finaldata\matcholder\with_pw\pword_matchedcon.txt'];
% [word,concrete,visual,hand,hear,valence,AoA,LEN1,FREQ1,OrthN,UN2_F,UN3_F]=textread(fname,'%s %f %f %f %f %f %f %f %f %f %f %f','headerlines',1);
% numpw=[concrete,visual,hand,hear,valence,AoA,LEN1,FREQ1,OrthN,UN2_F,UN3_F];
% 
% labels=[repmat({'visual'},50,1);repmat({'hear'},50,1);repmat({'hand'},50,1);repmat({'pseudoword'},50,1)];
% %word	1 concrete	2 visual	3 hand	4 hear	5 valence	6 AoA	7 LEN	8 FREQ	9 Orth	10 UN2_C	11UN3_C
% ylbl={'Concreteness', 'Visual', 'Hand-Action', 'Hearing',	'Valence (absolute difference from 5)',	'Age of Acquisition',	'Word Length',	'Word Frequency', 'Orthographic Neighbourhood',	'Bigram Frequency',	'Trigram Frequency'};
% %boxplot([numv(:,5);numhr(:,5);numhn(:,5);numna(:,5);numea(:,5)],labels)
% pmat=[];
% for ii=11
% 
% % anova1([numv(:,ii);numea(:,ii)],labels(1:100))
% p=kruskalwallis([numv(:,ii);numhr(:,ii);numhn(:,ii);numpw(:,ii)],labels)
% xlabel({'';'';'Category'})
% title(['p-value = ',num2str(p)])
% pmat(1)=kruskalwallis([numv(:,ii),numhr(:,ii)],[],'off');pmat(2)=kruskalwallis([numv(:,ii),numhn(:,ii)],[],'off');pmat(3)=kruskalwallis([numhr(:,ii),numhn(:,ii)],[],'off');
% pmat(4)=kruskalwallis([numv(:,ii),numpw(:,ii)],[],'off');pmat(5)=kruskalwallis([numhr(:,ii),numpw(:,ii)],[],'off');pmat(6)=kruskalwallis([numhn(:,ii),numpw(:,ii)],[],'off');
% ylabel(ylbl{ii})
% 
% end

% %% text read ca
% % 
% fname=['C:\Users\rf02\Documents\Rezvan\PhDproject\experimentdesign\nextexperiment\MEG\MEG_pilot\wordlist\word_concrete_matchedca50_2.txt'];
% [word,concrete,visual,hand,hear,valence,AoA,LEN1,FREQ1,OrthN,UN2_F,UN3_F]=textread(fname,'%s %f %f %f %f %f %f %f %f %f %f %f','headerlines',1);
% numc=[concrete,visual,hand,hear,valence,AoA,LEN1,FREQ1,OrthN,UN2_F,UN3_F];
% fname=['C:\Users\rf02\Documents\Rezvan\PhDproject\experimentdesign\nextexperiment\MEG\MEG_pilot\wordlist\word_eabs_matchedca50_2.txt'];
% [word,concrete,visual,hand,hear,valence,AoA,LEN1,FREQ1,OrthN,UN2_F,UN3_F]=textread(fname,'%s %f %f %f %f %f %f %f %f %f %f %f','headerlines',1);
% numea=[concrete,visual,hand,hear,valence,AoA,LEN1,FREQ1,OrthN,UN2_F,UN3_F];
% fname=['C:\Users\rf02\Documents\Rezvan\PhDproject\experimentdesign\nextexperiment\MEG\MEG_pilot\wordlist\word_nabs_matchedca50_2.txt'];
% [word,concrete,visual,hand,hear,valence,AoA,LEN1,FREQ1,OrthN,UN2_F,UN3_F]=textread(fname,'%s %f %f %f %f %f %f %f %f %f %f %f','headerlines',1);
% numna=[concrete,visual,hand,hear,valence,AoA,LEN1,FREQ1,OrthN,UN2_F,UN3_F];
% fname=['C:\Users\rf02\Documents\Rezvan\PhDproject\experimentdesign\nextexperiment\rating\finaldata\matcholder\with_pw\pword_matchedcon.txt'];
% [word,concrete,visual,hand,hear,valence,AoA,LEN1,FREQ1,OrthN,UN2_F,UN3_F]=textread(fname,'%s %f %f %f %f %f %f %f %f %f %f %f','headerlines',1);
% numpw=[concrete,visual,hand,hear,valence,AoA,LEN1,FREQ1,OrthN,UN2_F,UN3_F];
% numpw=numpw(1:40,:);
% 
% labels=[repmat({'concrete'},50,1);repmat({'neutral'},50,1);repmat({'emotional'},50,1)];%;repmat({'pseudoword'},45,1)];
% %word	1 concrete	2 visual	3 hand	4 hear	5 valence	6 AoA	7 LEN	8 FREQ	9 Orth	10 UN2_C	11UN3_C
% ylbl={'Concreteness', 'Visual', 'Hand-Action', 'Hearing',	'Valence (absolute difference from 5)',	'Age of Acquisition',	'Word Length',	'Word Frequency', 'Orthographic Neighbourhood',	'Bigram Frequency',	'Trigram Frequency'};
% %boxplot([numv(:,5);numhr(:,5);numhn(:,5);numna(:,5);numea(:,5)],labels)
% pmat=[];
% for ii=11
% % [h,p]=ttest2(numc(:,ii),numna(:,ii))
% % anova1([numv(:,ii);numea(:,ii)],labels(1:100))
% % anova1([numc(:,ii);numna(:,ii);numea(:,ii)],labels)
% p=kruskalwallis([numc(:,ii);numna(:,ii);numea(:,ii)],labels)%;numpw(:,ii)
% xlabel({'';'';'Category'})
% title(['p-value = ',num2str(p)])
% pmat(1)=kruskalwallis([numc(:,ii),numna(:,ii)],[],'off');pmat(2)=kruskalwallis([numc(:,ii),numea(:,ii)],[],'off');pmat(3)=kruskalwallis([numna(:,ii),numea(:,ii)],[],'off');
% % pmat(4)=kruskalwallis([numc(:,ii),numpw(:,ii)],[],'off');pmat(5)=kruskalwallis([numna(:,ii),numpw(:,ii)],[],'off');pmat(6)=kruskalwallis([numea(:,ii),numpw(:,ii)],[],'off');
% ylabel(ylbl{ii})
% end


% %% text read abs pw
% % 
% fname=['C:\Users\rf02\Documents\Rezvan\PhDproject\experimentdesign\nextexperiment\MEG\MEG_pilot\wordlist\pseudo_word_matchedabs.txt'];
% [word,concrete,visual,hand,hear,valence,AoA,LEN1,FREQ1,OrthN,UN2_F,UN3_F]=textread(fname,'%s %f %f %f %f %f %f %f %s %f %f %f','headerlines',1);
% numc=[concrete,visual,hand,hear,valence,AoA,LEN1,zeros(length(FREQ1),1),OrthN,UN2_F,UN3_F];
% fname=['C:\Users\rf02\Documents\Rezvan\PhDproject\experimentdesign\nextexperiment\MEG\MEG_pilot\wordlist\word_eabs_matchedpw.txt'];
% [word,concrete,visual,hand,hear,valence,AoA,LEN1,FREQ1,OrthN,UN2_F,UN3_F]=textread(fname,'%s %f %f %f %f %f %f %f %f %f %f %f','headerlines',1);
% numea=[concrete,visual,hand,hear,valence,AoA,LEN1,FREQ1,OrthN,UN2_F,UN3_F];
% fname=['C:\Users\rf02\Documents\Rezvan\PhDproject\experimentdesign\nextexperiment\MEG\MEG_pilot\wordlist\word_nabs_matchedpw.txt'];
% [word,concrete,visual,hand,hear,valence,AoA,LEN1,FREQ1,OrthN,UN2_F,UN3_F]=textread(fname,'%s %f %f %f %f %f %f %f %f %f %f %f','headerlines',1);
% numna=[concrete,visual,hand,hear,valence,AoA,LEN1,FREQ1,OrthN,UN2_F,UN3_F];
% fname=['C:\Users\rf02\Documents\Rezvan\PhDproject\experimentdesign\nextexperiment\rating\finaldata\matcholder\with_pw\pword_matchedcon.txt'];
% [word,concrete,visual,hand,hear,valence,AoA,LEN1,FREQ1,OrthN,UN2_F,UN3_F]=textread(fname,'%s %f %f %f %f %f %f %f %f %f %f %f','headerlines',1);
% numpw=[concrete,visual,hand,hear,valence,AoA,LEN1,FREQ1,OrthN,UN2_F,UN3_F];
% numpw=numpw(1:40,:);
% 
% labels=[repmat({'pword'},50,1);repmat({'neutral'},50,1);repmat({'emotional'},50,1)];%;repmat({'pseudoword'},45,1)];
% %word	1 concrete	2 visual	3 hand	4 hear	5 valence	6 AoA	7 LEN	8 FREQ	9 Orth	10 UN2_C	11UN3_C
% ylbl={'Concreteness', 'Visual', 'Hand-Action', 'Hearing',	'Valence (absolute difference from 5)',	'Age of Acquisition',	'Word Length',	'Word Frequency', 'Orthographic Neighbourhood',	'Bigram Frequency',	'Trigram Frequency'};
% %boxplot([numv(:,5);numhr(:,5);numhn(:,5);numna(:,5);numea(:,5)],labels)
% pmat=[];
% for ii=7
% % [h,p]=ttest2(numc(:,ii),numna(:,ii))
% % anova1([numv(:,ii);numea(:,ii)],labels(1:100))
% % anova1([numc(:,ii);numna(:,ii);numea(:,ii)],labels)
% p=kruskalwallis([numc(:,ii);numna(:,ii);numea(:,ii)],labels)%;numpw(:,ii)
% xlabel({'';'';'Category'})
% title(['p-value = ',num2str(p)])
% pmat(1)=kruskalwallis([numc(:,ii),numna(:,ii)],[],'off');pmat(2)=kruskalwallis([numc(:,ii),numea(:,ii)],[],'off');pmat(3)=kruskalwallis([numna(:,ii),numea(:,ii)],[],'off');
% % pmat(4)=kruskalwallis([numc(:,ii),numpw(:,ii)],[],'off');pmat(5)=kruskalwallis([numna(:,ii),numpw(:,ii)],[],'off');pmat(6)=kruskalwallis([numea(:,ii),numpw(:,ii)],[],'off');
% ylabel(ylbl{ii})
% end
% text read 3cat

% fname=['C:\Users\rf02\Documents\Rezvan\PhDproject\experimentdesign\nextexperiment\rating\finaldata\matcholder\3cat\word_vis_matched3.txt'];
fname=['C:\Users\rf02\Documents\Rezvan\PhDproject\experimentdesign\nextexperiment\MEG\MEG_pilot\wordlist\word_vis_matchedpw_afterMEG.txt'];
[word,concrete,visual,hand,hear,valence,AoA,LEN1,FREQ1,OrthN,UN2_F,UN3_F]=textread(fname,'%s %f %f %f %f %f %f %f %f %f %f %f','headerlines',1);
numv=[concrete,visual,hand,hear,valence,AoA,LEN1,FREQ1,OrthN,UN2_F,UN3_F];
% fname=['C:\Users\rf02\Documents\Rezvan\PhDproject\experimentdesign\nextexperiment\rating\finaldata\matcholder\3cat\word_hear_matched3.txt'];
fname=['C:\Users\rf02\Documents\Rezvan\PhDproject\experimentdesign\nextexperiment\MEG\MEG_pilot\wordlist\word_hear_matchedpw_afterMEG.txt'];
[word,concrete,visual,hand,hear,valence,AoA,LEN1,FREQ1,OrthN,UN2_F,UN3_F]=textread(fname,'%s %f %f %f %f %f %f %f %f %f %f %f','headerlines',1);
numhr=[concrete,visual,hand,hear,valence,AoA,LEN1,FREQ1,OrthN,UN2_F,UN3_F];
% fname=['C:\Users\rf02\Documents\Rezvan\PhDproject\experimentdesign\nextexperiment\rating\finaldata\matcholder\3cat\word_hand_matched3.txt'];
fname=['C:\Users\rf02\Documents\Rezvan\PhDproject\experimentdesign\nextexperiment\MEG\MEG_pilot\wordlist\word_hand_matchedpw_afterMEG.txt'];
[word,concrete,visual,hand,hear,valence,AoA,LEN1,FREQ1,OrthN,UN2_F,UN3_F]=textread(fname,'%s %f %f %f %f %f %f %f %f %f %f %f','headerlines',1);
numhn=[concrete,visual,hand,hear,valence,AoA,LEN1,FREQ1,OrthN,UN2_F,UN3_F];
fname=['C:\Users\rf02\Documents\Rezvan\PhDproject\experimentdesign\nextexperiment\MEG\MEG_pilot\wordlist\pseudo_word_matchedcon_afterMEG.txt'];
[word,concrete,visual,hand,hear,valence,AoA,LEN1,FREQ1,OrthN,UN2_F,UN3_F]=textread(fname,'%s %f %f %f %f %f %f %f %s %f %f %f','headerlines',1);
numpw=[concrete,visual,hand,hear,valence,AoA,LEN1,zeros(length(FREQ1),1),OrthN,UN2_F,UN3_F];

labels=[repmat({'visual'},49,1);repmat({'hear'},45,1);repmat({'hand'},50,1);repmat({'pseudoword'},50,1)];
%word	1 concrete	2 visual	3 hand	4 hear	5 valence	6 AoA	7 LEN	8 FREQ	9 Orth	10 UN2_C	11UN3_C
ylbl={'Concreteness', 'Visual', 'Hand-Action', 'Hearing',	'Valence (absolute difference from 5)',	'Age of Acquisition',	'Word Length',	'Word Frequency', 'Orthographic Neighbourhood',	'Bigram Frequency',	'Trigram Frequency'};
%boxplot([numv(:,5);numhr(:,5);numhn(:,5);numna(:,5);numea(:,5)],labels)
pmat=[];
for ii=1:11

% anova1([numv(:,ii);numea(:,ii)],labels(1:100))
p=kruskalwallis([numv(:,ii);numhr(:,ii);numhn(:,ii);numpw(:,ii)],labels)
xlabel({'';'';'Category'})
title(['p-value = ',num2str(p)])
% pmat(1)=kruskalwallis([numv(:,ii),numhr(:,ii)],[],'off');pmat(2)=kruskalwallis([numv(:,ii),numhn(:,ii)],[],'off');pmat(3)=kruskalwallis([numhr(:,ii),numhn(:,ii)],[],'off');
% pmat(4)=kruskalwallis([numv(:,ii),numpw(:,ii)],[],'off');pmat(5)=kruskalwallis([numhr(:,ii),numpw(:,ii)],[],'off');pmat(6)=kruskalwallis([numhn(:,ii),numpw(:,ii)],[],'off');
ylabel(ylbl{ii})

end

