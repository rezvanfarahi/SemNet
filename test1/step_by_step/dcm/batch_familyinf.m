% List of open inputs
nrun = X; % enter the number of runs here
jobfile = {'/home/rf02/rezvan/test1/step_by_step/dcm/batch_familyinf_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
end
spm('defaults', 'EEG');
spm_jobman('run', jobs, inputs{:});
