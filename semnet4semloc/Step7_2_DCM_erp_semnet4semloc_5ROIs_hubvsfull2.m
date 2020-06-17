clear
clc
close all
% addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/')
addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r6906/')

% addpath('\\cbsu\data\imaging\local\software\spm_cbu_svn\releases\spm12_latest\')

spm eeg
spm('defaults','eeg')
% dbstop if error
% clear

Model = struct([]);

%% Create models in advance
%order ATL IFG MTG AG vWFA
% Columns = FROM, rows = TO
%% hub ATL
%  For intrinsic (A) connections, {1} = forward; {2} = backward; {3} = lateral
Model(1).A{1} = [0 0 0 0 1
                1 0 0 0 0
                1 0 0 0 0
                1 0 0 0 0
                0 0 0 0 0];
Model(1).A{2} =[0 1 1 1 0
                0 0 0 0 0
                0 0 0 0 0
                0 0 0 0 0
                1 0 0 0 0];
Model(1).A{3} =[0 0 0 0 0
                0 0 0 0 0
                0 0 0 0 0
                0 0 0 0 0
                0 0 0 0 0];


% 1 For modulatory (B) connections, {n} = nth effect   , bidir-bidirs,
% input to WFA, no input self-modulation, spokes self-modulation
Model(1).B{1} =[1 1 1 1 1
                1 1 0 0 0
                1 0 1 0 0
                1 0 0 1 0
                1 0 0 0 0];
Model(1).C =    [0
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


%% hub MTG
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
Model(end+1) = Model(17);
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
Model(end+1) = Model(17);
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
Model(end+1) = Model(17);
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
Model(end+1) = Model(17);
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
Model(end+1)=Model(1);


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
Model(end+1)=Model(5);
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
Model(end+1)=Model(9);
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
Model(end+1)=Model(13);
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
length(Model)

%% Fit models

cwd = '/home/rf02/rezvan/semnet/semnet4semloc'

dosubs = [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19];%1 removed for LD
list_all =  {'./meg16_0030', 
            './meg16_0032', 
            './meg16_0034', 
            './meg16_0035', 
            './meg16_0042', 
            './meg16_0045', 
            './meg16_0052', 
            './meg16_0056',
            './meg16_0069', 
            './meg16_0070', 
            './meg16_0072', 
            './meg16_0073', 
            './meg16_0075', 
            './meg16_0078',
            './meg16_0082',
            './meg16_0086',
            './meg16_0097', 
            './meg16_0122', 
            './meg16_0125', 
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
this_model_dir='/imaging/rf02/Semnet/semnet4semloc/dcm/maxCTF/LD/oldreg_filtered_250ms_28models_5ROIs/maxCTF_ERP_dtr0/';%'./16models_200ms/semloc';
if ~exist(this_model_dir,'dir')
    mkdir(this_model_dir)
end
parfor s = 1:length(dosubs) %parfor
    ['subject ', num2str(s)]
    sub    = dosubs(s);
    
    DCMsub = DCMbase;
    fname_a='_Semnet_ConEmot_LD_oldreg_Signed_Evoked_5ROIs_meanflip_maxCTF_forDCM_avg';%'_SemLoc_Evoked_5ROIs_maxCTF_mean_forDCM_avg';%'_SemLoc_Evoked_5ROIs_ctfdip_avg_exttc_50verts_avg';%'_SemLoc_Evoked_5ROIs_meanCTF_50verts_aifg_forDCM_avg';%'_TypLex_Evoked_ATL_WFA_forDCM';
    path_in='/imaging/rf02/Semnet/semnet4semloc/dcm/maxCTF/';%'/imaging/rf02/TypLexMEG/dcm/5ROIs_hubvsfull/avg_50verts_dipole/';%'/imaging/rf02/TypLexMEG/dcm/5ROIs_whoshub/';
    fname_in=[path_in,list_all{sub}(3:end),fname_a,'_filtered_cnvrt'];
    fname_in
    %     this_dir='./5ROIs_whoshub/';
    %     fname_a='_SemLoc_Evoked_5ROIs_forDCM';
    %     fname_in=[list_all{s},fname_a,'_cnvrt'];
    %     spm_eeg_load(fname_in)
    S=[]; S.D = sprintf(fname_in)
    
    DCMsub.xY.Dfile = S.D;
    
    DCMsub = spm_dcm_erp_data(DCMsub,DCMsub.options.h);  % Either convert your LFP data into an SPM object or replace this function with something that returns DCM structure with relevant field
    
    DCMsub = spm_dcm_erp_dipfit(DCMsub);  % Either convert your LFP data into an SPM object or replace this function with something that returns DCM structure with relevant field
    
    for n=[2,4,39:40,6,8,41:42,10,12,43:44,14,16,45:46,18:2:36,37:38]%[22:2:36,37:38]%39:46%1:numel(Model)
        
        ['model ', num2str(n)]
        this_filename=sprintf([this_model_dir,'/DCM_erpf_ConEmot_LD_oldreg_5ROIs_avg_maxCTF_sub%d_mod%d.mat'],sub,n);
        if ~exist(this_filename,'file')
            DCM      = DCMsub;
            DCM.name = sprintf([this_model_dir,'/DCM_erpf_ConEmot_LD_oldreg_5ROIs_avg_maxCTF_sub%d_mod%d.mat'],sub,n);
            
            DCM.A = Model(n).A;
            DCM.B = Model(n).B;
            DCM.C = Model(n).C;
            %         DCM   = spm_dcm_erp(DCM);
            dcmcheck=0;
            while dcmcheck==0
                try
                    DCM   = spm_dcm_erp(DCM);
                    dcmcheck=1;
                end
            end
            
            
            %         LogEvd(s,n) = DCM.F;
            %         DCMname{s,n} = DCM.name;
        end
        
    end
end


%% Review results (have not cleaned)

% [dummy,maxi] = max(LogEvd);
%
% fprintf('Model %d is best',maxi);
%
% load(sprintf('DCM_erp_%d_%s',maxi,D.fname));
% % spm_dcm_erp_results_rik(DCM,'Coupling (C)');
% spm_dcm_erp_results_rik(DCM,'Coupling (A)');
% % spm_dcm_erp_results_rik(DCM,'Coupling (B)');
% spm_dcm_erp_results_rik(DCM,'Response');
%
% % spm_dcm_erp_results_rik(DCM,'ERPs (mode)');
% % spm_dcm_erp_results_rik(DCM,'ERPs (sources)');
%
%
% %    S=[]; S.D = sprintf('rtf_cbSPM_TFR_format_s%02d_retrieval_cr.mat',sub);
%     S=[]; S.D = sprintf('cbSPM_TFR_format_s%02d_retrieval_cr.mat',sub)
%
%     DCMsub.xY.Dfile = S.D;
%
%     DCMsub = spm_dcm_erp_dipfit(DCMsub, 1);
%
%     if ~isfield(DCMsub.xY,'source')  % (always true here)
%         DCMsub  = spm_dcm_ind_data(DCMsub);
%     end
%     DCMsub.options.gety = 0;
%     DCMsub.options.nograph  = 1;
%
%     for m=1:numel(model)
%
%         DCM      = DCMsub;
%         DCM.name = sprintf('DCM_ind_mod%d_sub%d_%s_%d-%d',m,sub,datafit,Twin(1),Twin(2));
%
%         DCM.A = model(m).A;
%         DCM.B = model(m).B;
%         DCM.C = model(m).C;
%
%         DCM   = spm_dcm_ind(DCM);   % saves automatically
%
%         LogEvd(s,m) = DCM.F;
%         DCMname{s,m} = DCM.name;
%     end
%
% end
%
