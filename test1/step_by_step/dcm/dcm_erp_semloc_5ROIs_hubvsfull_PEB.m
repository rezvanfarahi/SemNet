clear
clc
close all
addpath('/imaging/rf02/TypLexMEG/dcm/spm_12_latest/')%addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/')
% addpath('\\cbsu\data\imaging\local\software\spm_cbu_svn\releases\spm12_latest\')

spm eeg
spm('defaults','eeg')
% dbstop if error
% clear

Model = struct([]);

%% Create models in advance
%order ATL IFG MTG AG vWFA
% Columns = FROM, rows = TO
%% full model
Model(1).A{1} =[0 1 1 1 1
                1 0 1 1 1
                1 1 0 1 1
                1 1 1 0 1
                0 0 0 0 0];
Model(1).A{2} =[0 1 1 1 0
                1 0 1 1 0
                1 1 0 1 0
                1 1 1 0 0
                1 1 1 1 0];
Model(1).A{3} =[0 0 0 0 0
                0 0 0 0 0
                0 0 0 0 0
                0 0 0 0 0
                0 0 0 0 0];


% 1 For modulatory (B) connections, {n} = nth effect   , bidir-bidirs,
% input to WFA, no input self-modulation, spokes self-modulation
Model(1).B{1} =ones(5);
Model(1).C =    [0
    0
    0
    0
    1];

%% hub ATL
%  For intrinsic (A) connections, {1} = forward; {2} = backward; {3} = lateral
Model(end+1)=Model(end);
Model(end).A{1} =  [0 0 0 0 1
                    1 0 0 0 0
                    1 0 0 0 0
                    1 0 0 0 0
                    0 0 0 0 0];
Model(end).A{2} =  [0 1 1 1 0
                    0 0 0 0 0
                    0 0 0 0 0
                    0 0 0 0 0
                    1 0 0 0 0];
Model(end).A{3} =  [0 0 0 0 0
                    0 0 0 0 0
                    0 0 0 0 0
                    0 0 0 0 0
                    0 0 0 0 0];


% 1 For modulatory (B) connections, {n} = nth effect   , bidir-bidirs,
% input to WFA, no input self-modulation, spokes self-modulation
Model(end).B{1} =  [1 1 1 1 1
                    1 1 0 0 0
                    1 0 1 0 0
                    1 0 0 1 0
                    1 0 0 0 1];
Model(end).C =    [0
    0
    0
    0
    1];

% 2 For modulatory (B) connections, {n} = nth effect   , bidir-bidirs,
% input to WFA, input self-modulation, spoke self-modulation
Model(end+1)=Model(end);
Model(end).B{1} =  [1 1 1 1 1
                    1 1 0 0 0
                    1 0 1 0 0
                    1 0 0 1 0
                    1 0 0 0 1];


% 3 For modulatory (B) connections, {n} = nth effect   , bidir-bidirs,
% input to WFA, no input self-modulation, no spoke self-modulation
Model(end+1)=Model(end);
Model(end).B{1} =  [1 1 1 1 1
                    1 0 0 0 0
                    1 0 0 0 0
                    1 0 0 0 0
                    1 0 0 0 0];
% 4 For modulatory (B) connections, {n} = nth effect   , bidir-bidirs,
% input to WFA, input self-modulation, no spoke self-modulation
Model(end+1)=Model(end);
Model(end).B{1} =  [1 1 1 1 1
                    1 0 0 0 0
                    1 0 0 0 0
                    1 0 0 0 0
                    1 0 0 0 1];

%% hub IFG
% For intrinsic (A) connections, {1} = forward; {2} = backward; {3} = lateral
Model(end+1)=Model(end);
Model(end).A{1} =  [0 1 0 0 0
                    0 0 0 0 1
                    0 1 0 0 0
                    0 1 0 0 0
                    0 0 0 0 0];
Model(end).A{2} =  [0 0 0 0 0
                    1 0 1 1 0
                    0 0 0 0 0
                    0 0 0 0 0
                    0 1 0 0 0];
Model(end).A{3} =  [0 0 0 0 0
                    0 0 0 0 0
                    0 0 0 0 0
                    0 0 0 0 0
                    0 0 0 0 0];

% 5 For modulatory (B) connections, {n} = nth effect   , bidir-bidirs,
% input to WFA, no input self-modulation, spokes self-modulation
Model(end).B{1} =  [1 1 0 0 0
                    1 1 1 1 1
                    0 1 1 0 0
                    0 1 0 1 0
                    0 1 0 0 0];

Model(end).C =    [0
    0
    0
    0
    1];
% 6 For modulatory (B) connections, {n} = nth effect   , bidir-bidirs,
% input to WFA, input self-modulation, spokes self-modulation
Model(end+1)=Model(end);
Model(end).B{1} =  [1 1 0 0 0
                    1 1 1 1 1
                    0 1 1 0 0
                    0 1 0 1 0
                    0 1 0 0 1];


% 7 For modulatory (B) connections, {n} = nth effect   , bidir-bidirs,
% input to WFA, no input self-modulation, no spokes self-modulation
Model(end+1)=Model(end);
Model(end).B{1} =  [0 1 0 0 0
                    1 1 1 1 1
                    0 1 0 0 0
                    0 1 0 0 0
                    0 1 0 0 0];


% 8 For modulatory (B) connections, {n} = nth effect   ,bidir-bidirs,
% input to WFA, input self-modulation, no spokes self-modulation
Model(end+1)=Model(end);
Model(end).B{1} =  [0 1 0 0 0
                    1 1 1 1 1
                    0 1 0 0 0
                    0 1 0 0 0
                    0 1 0 0 1];


%% hub MFG
% For intrinsic (A) connections, {1} = forward; {2} = backward; {3} = lateral
Model(end+1)=Model(end);
Model(end).A{1} =  [0 0 1 0 0
                    0 0 1 0 0
                    0 0 0 0 1
                    0 0 1 0 0
                    0 0 0 0 0];
Model(end).A{2} = [0 0 0 0 0
                    0 0 0 0 0
                    1 1 0 1 0
                    0 0 0 0 0
                    0 0 1 0 0];
Model(end).A{3} = [0 0 0 0 0
                    0 0 0 0 0
                    0 0 0 0 0
                    0 0 0 0 0
                    0 0 0 0 0];

% 9 For modulatory (B) connections, {n} = nth effect   , bidir-bidirs,
% input to WFA, no input self-modulation, spokes self-modulation
Model(end).B{1} =  [1 0 1 0 0
                    0 1 1 0 0
                    1 1 1 1 1
                    0 0 1 1 0
                    0 0 1 0 0];

Model(end).C =    [0
    0
    0
    0
    1];
% 10 For modulatory (B) connections, {n} = nth effect   , bidir-bidirs,
% input to WFA, input self-modulation, spokes self-modulation
Model(end+1)=Model(end);
Model(end).B{1} =  [1 0 1 0 0
                    0 1 1 0 0
                    1 1 1 1 1
                    0 0 1 1 0
                    0 0 1 0 1];


% 11 For modulatory (B) connections, {n} = nth effect   , bidir-bidirs,
% input to WFA, no input self-modulation, no spokes self-modulation
Model(end+1)=Model(end);
Model(end).B{1} =  [0 0 1 0 0
                    0 0 1 0 0
                    1 1 1 1 1
                    0 0 1 0 0
                    0 0 1 0 0];


% 12 For modulatory (B) connections, {n} = nth effect   , bidir-bidirs,
% input to WFA,  input self-modulation, no spokes self-modulation
Model(end+1)=Model(end);
Model(end).B{1} =  [0 0 1 0 0
                    0 0 1 0 0
                    1 1 1 1 1
                    0 0 1 0 0
                    0 0 1 0 1];



%% hub AG
% For intrinsic (A) connections, {1} = forward; {2} = backward; {3} = lateral
Model(end+1)=Model(end);
Model(end).A{1} =  [0 0 0 1 0
                    0 0 0 1 0
                    0 0 0 1 0
                    0 0 0 0 1
                    0 0 0 0 0];
Model(end).A{2} =  [0 0 0 0 0
                    0 0 0 0 0
                    0 0 0 0 0
                    1 1 1 0 0
                    0 0 0 1 0];
Model(end).A{3} =  [0 0 0 0 0
    0 0 0 0 0
    0 0 0 0 0
    0 0 0 0 0
    0 0 0 0 0];

% 13 For modulatory (B) connections, {n} = nth effect   , bidir-bidirs,
% input to WFA, no input self-modulation, spokes self-modulation
Model(end).B{1} =  [1 0 0 1 0
                    0 1 0 1 0
                    0 0 1 1 0
                    1 1 1 1 1
                    0 0 0 1 0];

Model(end).C =    [0
    0
    0
    0
    1];
% 14 For modulatory (B) connections, {n} = nth effect   , bidir-bidirs,
% input to WFA, input self-modulation, spokes self-modulation
Model(end+1)=Model(end);
Model(end).B{1} =  [1 0 0 1 0
                    0 1 0 1 0
                    0 0 1 1 0
                    1 1 1 1 1
                    0 0 0 1 1];


% 15 For modulatory (B) connections, {n} = nth effect   ,  bidir-bidirs,
% input to WFA, no input self-modulation, no spokes self-modulation
Model(end+1)=Model(end);
Model(end).B{1} =  [0 0 0 1 0
                    0 0 0 1 0
                    0 0 0 1 0
                    1 1 1 1 1
                    0 0 0 1 0];

% 16 For modulatory (B) connections, {n} = nth effect   , bidir-bidirs,
% input to WFA,  input self-modulation, no spokes self-modulation
Model(end+1)=Model(end);
Model(end).B{1} =  [0 0 0 1 0
                    0 0 0 1 0
                    0 0 0 1 0
                    1 1 1 1 1
                    0 0 0 1 1];


%% convergence zones
% For intrinsic (A) connections, {1} = forward; {2} = backward; {3} = lateral
Model(end+1) = Model(1);
Model(end).A{1} =  [0 0 0 0 1
                    0 0 0 0 1
                    0 0 0 0 1
                    0 0 0 0 1
                    0 0 0 0 0];
Model(end).A{2} =  [0 0 0 0 0
                    0 0 0 0 0
                    0 0 0 0 0
                    0 0 0 0 0
                    1 1 1 1 0];
Model(end).A{3} = [0 0 0 0 0
    0 0 0 0 0
    0 0 0 0 0
    0 0 0 0 0
    0 0 0 0 0];
%% ALL
% 17 For modulatory (B) connections, {n} = nth effect   , bidirs,
% input to WFA, no input self-modulation, spokes self-modulation
Model(end).B{1} =  [1 0 0 0 1
                    0 1 0 0 1
                    0 0 1 0 1
                    0 0 0 1 1
                    1 1 1 1 0];

Model(end).C =    [0
    0
    0
    0
    1];
% 18 For modulatory (B) connections, {n} = nth effect   , bidirs,
% input to WFA, input self-modulation, spokes self-modulation


Model(end+1) = Model(end);
Model(end).B{1} =  [1 0 0 0 1
                    0 1 0 0 1
                    0 0 1 0 1
                    0 0 0 1 1
                    1 1 1 1 1];
Model(end).C =    [0
    0
    0
    0
    1];
% 19 For modulatory (B) connections, {n} = nth effect   , bidirs,
% input to WFA, no input self-modulation, no spokes self-modulation


Model(end+1) = Model(end);
Model(end).B{1} =  [0 0 0 0 1
                    0 0 0 0 1
                    0 0 0 0 1
                    0 0 0 0 1
                    1 1 1 1 0];
Model(end).C =    [0
    0
    0
    0
    1];
% 20 For modulatory (B) connections, {n} = nth effect   , bidirs,
% input to WFA,  input self-modulation, no spokes self-modulation


Model(end+1) = Model(end);
Model(end).B{1} =  [0 0 0 0 1
                    0 0 0 0 1
                    0 0 0 0 1
                    0 0 0 0 1
                    1 1 1 1 1];
Model(end).C =    [0
    0
    0
    0
    1];

%% ATL
% 21 For modulatory (B) connections, {n} = nth effect   , bidirs,
% input to WFA, no input self-modulation, spokes self-modulation
Model(end+1) = Model(17+1);
Model(end).B{1} =  [1 0 0 0 1
                    0 0 0 0 0
                    0 0 0 0 0
                    0 0 0 0 0
                    1 0 0 0 0];

Model(end).C =    [0
    0
    0
    0
    1];
% 22 For modulatory (B) connections, {n} = nth effect   , bidirs,
% input to WFA, input self-modulation, spokes self-modulation


Model(end+1) = Model(end);
Model(end).B{1} =  [1 0 0 0 1
                    0 0 0 0 0
                    0 0 0 0 0
                    0 0 0 0 0
                    1 0 0 0 1];
Model(end).C =    [0
    0
    0
    0
    1];
% 23 For modulatory (B) connections, {n} = nth effect   , bidirs,
% input to WFA, no input self-modulation, no spokes self-modulation


Model(end+1) = Model(end);
Model(end).B{1} =  [0 0 0 0 1
                    0 0 0 0 0
                    0 0 0 0 0
                    0 0 0 0 0
                    1 0 0 0 0];
Model(end).C =    [0
    0
    0
    0
    1];
% 24 For modulatory (B) connections, {n} = nth effect   , bidirs,
% input to WFA,  input self-modulation, no spokes self-modulation


Model(end+1) = Model(end);
Model(end).B{1} =  [0 0 0 0 1
                    0 0 0 0 0
                    0 0 0 0 0
                    0 0 0 0 0
                    1 0 0 0 1];
Model(end).C =    [0
    0
    0
    0
    1];
%% IFG
% 25 For modulatory (B) connections, {n} = nth effect   , bidirs,
% input to WFA, no input self-modulation, spokes self-modulation
Model(end+1) = Model(17+1);
Model(end).B{1} =  [0 0 0 0 0
                    0 1 0 0 1
                    0 0 0 0 0
                    0 0 0 0 0
                    0 1 0 0 0];

Model(end).C =    [0
    0
    0
    0
    1];
% 26 For modulatory (B) connections, {n} = nth effect   , bidirs,
% input to WFA, input self-modulation, spokes self-modulation


Model(end+1) = Model(end);
Model(end).B{1} =  [0 0 0 0 0
                    0 1 0 0 1
                    0 0 0 0 0
                    0 0 0 0 0
                    0 1 0 0 1];
Model(end).C =    [0
    0
    0
    0
    1];
% 27 For modulatory (B) connections, {n} = nth effect   , bidirs,
% input to WFA, no input self-modulation, no spokes self-modulation


Model(end+1) = Model(end);
Model(end).B{1} =  [0 0 0 0 0
                    0 0 0 0 1
                    0 0 0 0 0
                    0 0 0 0 0
                    0 1 0 0 0];
Model(end).C =    [0
    0
    0
    0
    1];
% 28 For modulatory (B) connections, {n} = nth effect   , bidirs,
% input to WFA,  input self-modulation, no spokes self-modulation


Model(end+1) = Model(end);
Model(end).B{1} =  [0 0 0 0 0
                    0 0 0 0 1
                    0 0 0 0 0
                    0 0 0 0 0
                    0 1 0 0 1];
Model(end).C =    [0
    0
    0
    0
    1];
%% MTG
% 29 For modulatory (B) connections, {n} = nth effect   , bidirs,
% input to WFA, no input self-modulation, spokes self-modulation
Model(end+1) = Model(17+1);
Model(end).B{1} =  [0 0 0 0 0
                    0 0 0 0 0
                    0 0 1 0 1
                    0 0 0 0 0
                    0 0 1 0 0];

Model(end).C =    [0
    0
    0
    0
    1];
% 30 For modulatory (B) connections, {n} = nth effect   , bidirs,
% input to WFA, input self-modulation, spokes self-modulation


Model(end+1) = Model(end);
Model(end).B{1} =  [0 0 0 0 0
                    0 0 0 0 0
                    0 0 1 0 1
                    0 0 0 0 0
                    0 0 1 0 1];
Model(end).C =    [0
    0
    0
    0
    1];
% 31 For modulatory (B) connections, {n} = nth effect   , bidirs,
% input to WFA, no input self-modulation, no spokes self-modulation


Model(end+1) = Model(end);
Model(end).B{1} =  [0 0 0 0 0
                    0 0 0 0 0
                    0 0 0 0 1
                    0 0 0 0 0
                    0 0 1 0 0];
Model(end).C =    [0
    0
    0
    0
    1];
% 32 For modulatory (B) connections, {n} = nth effect   , bidirs,
% input to WFA,  input self-modulation, no spokes self-modulation


Model(end+1) = Model(end);
Model(end).B{1} =  [0 0 0 0 0
                    0 0 0 0 0
                    0 0 0 0 1
                    0 0 0 0 0
                    0 0 1 0 1];
Model(end).C =    [0
    0
    0
    0
    1];
%% AG
% 33 For modulatory (B) connections, {n} = nth effect   , bidirs,
% input to WFA, no input self-modulation, spokes self-modulation
Model(end+1) = Model(17+1);
Model(end).B{1} =  [0 0 0 0 0
                    0 0 0 0 0
                    0 0 0 0 0
                    0 0 0 1 1
                    0 0 0 1 0];

Model(end).C =    [0
    0
    0
    0
    1];
% 34 For modulatory (B) connections, {n} = nth effect   , bidirs,
% input to WFA, input self-modulation, spokes self-modulation


Model(end+1) = Model(end);
Model(end).B{1} =  [0 0 0 0 0
                    0 0 0 0 0
                    0 0 0 0 0
                    0 0 0 1 1
                    0 0 0 1 1];
Model(end).C =    [0
    0
    0
    0
    1];
% 35 For modulatory (B) connections, {n} = nth effect   , bidirs,
% input to WFA, no input self-modulation, no spokes self-modulation


Model(end+1) = Model(end);
Model(end).B{1} =  [0 0 0 0 0
                    0 0 0 0 0
                    0 0 0 0 0
                    0 0 0 0 1
                    0 0 0 1 0];
Model(end).C =    [0
    0
    0
    0
    1];
% 36 For modulatory (B) connections, {n} = nth effect   , bidirs,
% input to WFA,  input self-modulation, no spokes self-modulation


Model(end+1) = Model(end);
Model(end).B{1} =  [0 0 0 0 0
                    0 0 0 0 0
                    0 0 0 0 0
                    0 0 0 0 1
                    0 0 0 1 1];
Model(end).C =    [0
    0
    0
    0
    1];
% 37 Null (baseline) model 1 - no modulation by Fac input self-modulation

Model(end+1) = Model(end);
Model(end).A{1} = [0 0 0 0 0
    0 0 0 0 0
    0 0 0 0 0
    0 0 0 0 0
    0 0 0 0 0];
Model(end).A{2} = [0 0 0 0 0
    0 0 0 0 0
    0 0 0 0 0
    0 0 0 0 0
    0 0 0 0 0];
Model(end).A{3} = [0 0 0 0 0
    0 0 0 0 0
    0 0 0 0 0
    0 0 0 0 0
    0 0 0 0 0];
Model(end).B{1} = [0 0 0 0 0
    0 0 0 0 0
    0 0 0 0 0
    0 0 0 0 0
    0 0 0 0 1];
% 38 Null (baseline) model 1 - no modulation by Fac , no input self-modulation

Model(end+1) = Model(end);
Model(end).B{1} = zeros(5);
%% hub ATL-continued
%  For intrinsic (A) connections, {1} = forward; {2} = backward; {3} = lateral
Model(end+1)=Model(1+1);


% 39 For modulatory (B) connections, {n} = nth effect   , bidir-bidirs,
% input to WFA, input self-modulation, spokes self-modulation, no hub self
% modulation
Model(end).B{1} =  [0 1 1 1 1
                    1 1 0 0 0
                    1 0 1 0 0
                    1 0 0 1 0
                    1 0 0 0 1];
Model(end).C =    [0
    0
    0
    0
    1];

% 40 For modulatory (B) connections, {n} = nth effect   , bidir-bidirs,
% input to WFA, input self-modulation, no spokes self-modulation, no hub self
% modulation
Model(end+1) = Model(end);
Model(end).B{1} =  [0 1 1 1 1
                    1 0 0 0 0
                    1 0 0 0 0
                    1 0 0 0 0
                    1 0 0 0 1];

%% hub IFG continued
% For intrinsic (A) connections, {1} = forward; {2} = backward; {3} = lateral
Model(end+1)=Model(5+1);
Model(end).A{1} =  [0 1 0 0 0
                    0 0 0 0 1
                    0 1 0 0 0
                    0 1 0 0 0
                    0 0 0 0 0];
Model(end).A{2} =  [0 0 0 0 0
                    1 0 1 1 0
                    0 0 0 0 0
                    0 0 0 0 0
                    0 1 0 0 0];
Model(end).A{3} =  [0 0 0 0 0
                    0 0 0 0 0
                    0 0 0 0 0
                    0 0 0 0 0
                    0 0 0 0 0];

% 41 For modulatory (B) connections, {n} = nth effect   , bidir-bidirs,
% input to WFA, input self-modulation, spokes self-modulation, no hub self
% modulation
Model(end).B{1} =  [1 1 0 0 0
                    1 0 1 1 1
                    0 1 1 0 0
                    0 1 0 1 0
                    0 1 0 0 1];

Model(end).C =    [0
    0
    0
    0
    1];
% 42 For modulatory (B) connections, {n} = nth effect   , bidir-bidirs,
% input to WFA, input self-modulation, no spokes self-modulation, no hub self
% modulation
Model(end+1)=Model(end);
Model(end).B{1} =  [0 1 0 0 0
                    1 0 1 1 1
                    0 1 0 0 0
                    0 1 0 0 0
                    0 1 0 0 1];

%% hub MFG continued
% For intrinsic (A) connections, {1} = forward; {2} = backward; {3} = lateral
Model(end+1)=Model(9+1);
Model(end).A{1} =  [0 0 1 0 0
                    0 0 1 0 0
                    0 0 0 0 1
                    0 0 1 0 0
                    0 0 0 0 0];
Model(end).A{2} =  [0 0 0 0 0
                    0 0 0 0 0
                    1 1 0 1 0
                    0 0 0 0 0
                    0 0 1 0 0];
Model(end).A{3} =  [0 0 0 0 0
                    0 0 0 0 0
                    0 0 0 0 0
                    0 0 0 0 0
                    0 0 0 0 0];

% 43 For modulatory (B) connections, {n} = nth effect   , bidir-bidirs,
% input to WFA, input self-modulation, spokes self-modulation, no hub self
% modulation
Model(end).B{1} =  [1 0 1 0 0
                    0 1 1 0 0
                    1 1 0 1 1
                    0 0 1 1 0
                    0 0 1 0 1];

Model(end).C =    [0
    0
    0
    0
    1];
% 44 For modulatory (B) connections, {n} = nth effect   , bidir-bidirs,
% input to WFA, input self-modulation, no spokes self-modulation, no hub self
% modulation
Model(end+1)=Model(end);
Model(end).B{1} =  [0 0 1 0 0
                    0 0 1 0 0
                    1 1 0 1 1
                    0 0 1 0 0
                    0 0 1 0 1];

%% hub AG continued
% For intrinsic (A) connections, {1} = forward; {2} = backward; {3} = lateral
Model(end+1)=Model(13+1);
Model(end).A{1} =  [0 0 0 1 0
                    0 0 0 1 0
                    0 0 0 1 0
                    0 0 0 0 1
                    0 0 0 0 0];
Model(end).A{2} =  [0 0 0 0 0
                    0 0 0 0 0
                    0 0 0 0 0
                    1 1 1 0 0
                    0 0 0 1 0];
Model(end).A{3} =  [0 0 0 0 0
    0 0 0 0 0
    0 0 0 0 0
    0 0 0 0 0
    0 0 0 0 0];

% 45 For modulatory (B) connections, {n} = nth effect   , bidir-bidirs,
% input to WFA, input self-modulation, spokes self-modulation, no hub self
% modulation
Model(end).B{1} =  [1 0 0 1 0
                    0 1 0 1 0
                    0 0 1 1 0
                    1 1 1 0 1
                    0 0 0 1 1];

Model(end).C =    [0
    0
    0
    0
    1];
% 46 For modulatory (B) connections, {n} = nth effect   , bidir-bidirs,
% input to WFA, input self-modulation, no spokes self-modulation, no hub self
% modulation
Model(end+1)=Model(end);
Model(end).B{1} =  [0 0 0 1 0
                    0 0 0 1 0
                    0 0 0 1 0
                    1 1 1 0 1
                    0 0 0 1 1];
             

%% Fit models

cwd = '/home/rf02/rezvan/test1/dcm/5ROIs_whoshub'

dosubs = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17];
list_all = {'./meg10_0378',
    './meg10_0390',
    './meg11_0026',
    './meg11_0050',
    './meg11_0052',
    './meg11_0069',
    './meg11_0086',
    './meg11_0091',
    './meg11_0096',
    './meg11_0101',
    './meg11_0102',
    './meg11_0112',
    './meg11_0104',
    './meg11_0118',
    './meg11_0131',
    './meg11_0144',
    './meg11_0147',
    };
% Set up DCM structure

DCMbase = [];
DCMbase.xY.modality = 'LFP';
DCMbase.options.analysis = 'ERP';
DCMbase.options.model    = 'ERP';
DCMbase.options.spatial =  'LFP';
DCMbase.options.trials  = [1 2];  % [3 4] for CROI ([1 2] ROI)
DCMbase.options.Tdcm(1) = 0;      % start of peri-stimulus time to be modelled
DCMbase.options.Tdcm(2) = 250;   % end of peri-stimulus time to be modelled
DCMbase.options.Nmodes  = 8;      % nr of modes for data selection
DCMbase.options.h       = 0;      % nr of DCT components for detrending [increase to 4????]
DCMbase.options.onset   = 60;     % selection of onset (prior mean = 60ms for onset at 0)
DCMbase.options.D       = 1;      % downsampling
DCMbase.options.han     = 0;      % Hanning

DCMbase.options.lock     = 0;      % (if want modulations to be in same direction for all connections)
DCMbase.options.location = 0;      % (Not relevant; only for ECD)
DCMbase.options.symmetry = 0;      % (Not relevant; only for ECD)

DCMbase.options.nograph  = 0;

% Region names (Location, Lpos, irrelevant for LFP)
DCMbase.Sname = {'lATL','lIFG','lMTG','lAG','lWFA'};
DCMbase.Lpos = [];

DCMbase.xU.X = [1; -1]; % W-NW
DCMbase.xU.name = {'C-A'};


% LogEvd=[]; DCMname={};
this_model_dir='/imaging/rf02/TypLexMEG/dcm/5ROIs_hubvsfull/filtered_250ms_28models_5ROIs/semloc/maxCTF_ERP_dtr0_PEB/fitbmr/';%'/imaging/rf02/TypLexMEG/dcm/5ROIs_hubvsfull/250ms_38models_5ROIs/semloc0';%'./16models_200ms/semloc';
if ~exist(this_model_dir,'dir')
    mkdir(this_model_dir)
end
%% create GCM
GCM={};
for s = 1:length(dosubs) %parfor
    s
    sub    = dosubs(s);
    
    DCMsub = DCMbase;
    fname_a='_SemLoc_Evoked_5ROIs_maxCTF_mean_forDCM_avg';%'_SemLoc_Evoked_5ROIs_maxCTF_forDCM_avg';%'_SemLoc_Evoked_5ROIs_meanCTF_forDCM';%'_TypLex_Evoked_ATL_WFA_forDCM';
    path_in='/imaging/rf02/TypLexMEG/dcm/maxCTF/';%'/imaging/rf02/TypLexMEG/dcm/5ROIs_hubvsfull/';%'/imaging/rf02/TypLexMEG/dcm/5ROIs_whoshub/';
    fname_in=[path_in,list_all{sub}(3:end),fname_a,'_filtered_cnvrt'];
    %     this_dir='./5ROIs_whoshub/';
    %     fname_a='_SemLoc_Evoked_5ROIs_forDCM';
    %     fname_in=[list_all{s},fname_a,'_cnvrt'];
    %     spm_eeg_load(fname_in)
    S=[]; S.D = sprintf(fname_in)
    
    DCMsub.xY.Dfile = S.D;
    
    DCMsub = spm_dcm_erp_data(DCMsub,DCMsub.options.h);  % Either convert your LFP data into an SPM object or replace this function with something that returns DCM structure with relevant field
    
    DCMsub = spm_dcm_erp_dipfit(DCMsub);  % Either convert your LFP data into an SPM object or replace this function with something that returns DCM structure with relevant field
    ncnt=0;
    for n=[0,2,4,39:40,6,8,41:42,10,12,43:44,14,16,45:46,18:2:36,37:38]+1%39:46%1:numel(Model)
        ncnt=ncnt+1;
        n
        this_filename=sprintf([this_model_dir,'/DCMPEB_erp_SemLoc_maxCTF_5ROIs_sub%d_mod%d.mat'],sub,n);
            DCM      = DCMsub;
            DCM.name = sprintf([this_model_dir,'/DCMPEB_erp_SemLoc_maxCTF_5ROIs_sub%d_mod%d.mat'],sub,n);
            
            DCM.A = Model(n).A;
            DCM.B = Model(n).B;
            DCM.C = Model(n).C;
            GCM{s,ncnt}=DCM;
            %         DCM   = spm_dcm_erp(DCM);
%             dcmcheck=0;
%             while dcmcheck==0
%                 try
%                     DCM   = spm_dcm_erp(DCM);
%                     dcmcheck=1;
%                 end
%             end
            
            
            %         LogEvd(s,n) = DCM.F;
            %         DCMname{s,n} = DCM.name;
        
    end
end
%% DCM fit and BMR
%%% fit most complex first, then bmr- start
GCM1=GCM;
outpath=[this_model_dir,'/GCM_fwdbwd_fitbmr_modeldefs.mat'];%'/imaging/rf02/TypLexMEG/dcm/5ROIs_hubvsfull/250ms_38models_5ROIs/semloc0/GCM_fwdbwd_fitbmr_modeldefs.mat';
save(outpath,'GCM')
GCM(:,1) = spm_dcm_fit(GCM(:,1)); 
GCM2=GCM;
outpath=[this_model_dir,'/GCM_fwdbwd_fitbmr_fullestimated.mat'];%'/imaging/rf02/TypLexMEG/dcm/5ROIs_hubvsfull/250ms_38models_5ROIs/semloc0/GCM_fwdbwd_fitbmr_fullestimated.mat';
save(outpath,'GCM')
[RCM,BMC,BMA] = spm_dcm_bmr(GCM2);
GCM3={};
for ii = 1:length(RCM)
    for jj = 1:length(RCM{1})
        GCM3{ii,jj}=RCM{ii,1}{1,jj};
    end
end
GCM=GCM3;
outpath=[this_model_dir,'/GCM_fwdbwd_fitbmr_RCM_BMR.mat'];%'/imaging/rf02/TypLexMEG/dcm/5ROIs_hubvsfull/250ms_38models_5ROIs/semloc0/GCM_fwdbwd_fitbmr_RCM_BMR.mat';
save(outpath,'GCM')        
%%% fit most complex first, then bmr- finish

%%% fit_peb, draw subjects from local minima- start
% GCM1=GCM;
% outpath=[this_model_dir,'/GCM_fwdbwd_pebfit_modeldefs.mat'];%'/imaging/rf02/TypLexMEG/dcm/5ROIs_hubvsfull/250ms_38models_5ROIs/semloc0/GCM_fwdbwd_pebfit_modeldefs.mat';
% save(outpath,'GCM')
% GCM = spm_dcm_peb_fit(GCM);  
% GCM2=GCM;
% outpath=[this_model_dir,'/GCM_fwdbwd_pebfit_fullestimated.mat'];%outpath='/imaging/rf02/TypLexMEG/dcm/5ROIs_hubvsfull/250ms_38models_5ROIs/semloc0/GCM_fwdbwd_pebfit_fullestimated.mat';
% save(outpath,'GCM')
% GCM3={};
% for ii = 1:length(GCM)
%     for jj = 1:length(GCM{1})
%         GCM3{ii,jj}=GCM{ii,1}{1,jj};
%     end
% end
% GCM=GCM3;
% outpath=[this_model_dir,'/GCM_fwdbwd_pebfit_RCM_BMR.mat'];%'/imaging/rf02/TypLexMEG/dcm/5ROIs_hubvsfull/250ms_38models_5ROIs/semloc0/GCM3_fwdbwd_pebfit_RCM_BMR.mat';
% save(outpath,'GCM')
%%% fit_peb, draw subjects from local minima- finish

% F = zeros(size(GCM,1),size(GCM,2)-1);
% for s = 1:size(GCM,1)
%     for m = 2:size(GCM,2);
%         F(s,m-1) = GCM{s,m}.F;
%     end
% end
% save('log_evidence_matrix_GCM3_fwdbwd_pebfit.mat', 'F');

%% running PEB 
M=struct();
M.alpha=1;
M.beta=16;%
M.hE    = 0;
M.hC    = 1/16;
M.Q     = 'single';
M.X = ones(17,1);
% Choose field
field = {'A','B'};
% Estimate model
PEB     = spm_dcm_peb(GCM,M,field);
outpath=[this_model_dir,'/PEB_fwdbwd_fitbmr_RCM_BMR.mat'];
load(outpath,'PEB')
BMA2 = spm_dcm_peb_bmc(PEB(1),GCM(1,:));%spm_dcm_peb_bmc(PEB(1));

spm_dcm_peb_review(PEB(1),GCM(1,:))
