clear S
addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest')
DIR = '/imaging/ta02/ActionWords/LD/aamod_firstlevel_contrasts_00001';
[junk, uname] = system('whoami'); uname = uname(1:end-1);
S.outdir = ['/home/' uname '/rezvan/test1'];

% Order: Easy_WORD (4), Easy_NONWORD (3), Hard_WORD (8), Hard_NONWORD (7)
d = dir(fullfile(DIR,'CBU*'));
males = [2 5 6 7 10 12 17]; im = 0;
females = [1 3 4 8 9 11 13 14 15 16 18 19]; ifm = 0;
for i = 1:numel(d)
    fcons = spm_select('FPListRec',fullfile(DIR,d(i).name),'con_000.\.img');   
    if ~isempty(find(i==males,1))
        im = im + 1;
        S.imgfiles{1}{im} = fcons([7],:);
    elseif ~isempty(find(i==females,1))
        ifm = ifm + 1;
        S.imgfiles{2}{ifm} = fcons([7],:);
    end
end

S.contrasts{1}.name = 'Task G1';
S.contrasts{1}.type = 'T';
S.contrasts{1}.c = [1 0];

S.contrasts{2}.name = 'Task G2';
S.contrasts{2}.type = 'T';
S.contrasts{2}.c = [0 1];

S.contrasts{3}.name = 'Task';
S.contrasts{3}.type = 'T';
S.contrasts{3}.c = [1 1];

S.contrasts{4}.name = 'Task G1>G2';
S.contrasts{4}.type = 'T';
S.contrasts{4}.c = [1 -1];


% contrasts   - cell array of contrast structures, with fields c
%                  (matrix), type ('F' or 'T') and name (optional)
spm('Defaults','FMRI');
batch_spm_anova(S);



