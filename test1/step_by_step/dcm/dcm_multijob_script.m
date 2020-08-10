% This script is intended to provide a simple illustration of how multiple
% jobs can be submitted to the CBU queueing system using Matlab. In this
% case, a "job" is defined as any task that you want a single worker node
% (equivalent to a single CPU core on our cluster) to work on. 

% How you divide your analysis into different jobs will depend on several 
% factors, including any dependencies that may exist between different
% elements of your analysis (e.g. does the data need to be spatially
% aligned before it can be normalised, etc). In a typical fMRI analysis,
% one obvious way to divide an analysis into jobs would be to treat the
% analysis of each subject's data as a single job.

% You could also sub-divide the analysis into even smaller chunks - e.g.
% by splitting the pre-processing steps for each subject into separate
% jobs. In this case however, there would be dependencies between the jobs,
% and you would need to make the scheduler aware of those dependencies. If
% you use the cbu_qsub wrapper function (see code below), you'll create a
% structure array containing details of all your jobs. This jobs structure
% contains a 'depends_on' field where, for each job, you can specify which
% other jobs need to run first.
% 

% Once you have decided how to divide your analysis, you need to create a
% job-specific scripts (or functions) to run each individual job. These
% can be written to take input arguments, so one script could 
% handle many differnt jobs - e.g. my_preprocessing_script('subject-1'),
% my_preprocessing_script('subject-2'), etc (or even 
% my_preprocessing_script('realign','subject-1'), etc if you wanted to
% split jobs up further). The current script calls 2 example single job
% scripts - example_singlejob_preprocessing, and
% example_singlejob_firstlevel. Both of these scripts have been written to
% take a single input argument, subject id.


%% ========================================================================
% specify various parameters

% list of subjects to analyse:
subjects=dir('/imaging/russell/test_DCS/CBU*');
subjects={subjects(:).name};
subjects=subjects(1:3);
subjects = {'meg10_0378',
    'meg10_0390',
    'meg11_0026',
    'meg11_0050',
    'meg11_0052',
    'meg11_0069',
    'meg11_0086',
    'meg11_0091',
    'meg11_0096',
    'meg11_0101',
    'meg11_0102',
    'meg11_0112',
    'meg11_0104',
    'meg11_0118',
    'meg11_0131',
    'meg11_0144',
    'meg11_0147',
    };
% list of directories that should be added to the matlab path on worker
% nodes. Worker nodes are fresh instances of matlab, and will start with 
% your default matlab path. If your job scripts are not on your default 
% path, include their directories below:
workerpath={'/home/rf02/rezvan/test1/step_by_step/dcm/'};

% add contents of workerpath to the current matlab path
addpath(workerpath{:},'-BEGIN');

% list of tasks to run for each subject - for each subject, run the
% pre-processing, then the first level model.
tasks={'single_dcm_erp_typlex', ...
    'single_dcm_erp_semloc'};

% directory where analysis results will be stored
anadir='/imaging/rf02/TypLexMEG/dcm/2ROIs_whoshub';



%% ========================================================================
% Create a job array containing a list of all jobs the jobs you want to run
% in the format expected by the cbu_qsub function. This is a struct array 
% with fields:
%   .task - a function handle for the script that will run this task (e.g.
%       @/hpc-software/matlab/cbu/example_singlejob_firstlevel)
%   .n_return_values - integer, number of values returned by this task
%       script
%   .input_args - cell array - one element for each input argument to be
%       passed to the task script
%   .depends_on - vector of numbers representing the jobs that must
%       complete before running the current job. Jobs are represented by
%       their index number in the job array.

% In this example, we are running the same script for multiple datasets, so
% we can build the job array using 2 loops:


njob=size(tasks,2)*size(subjects,1); % we are creating 2 jobs for each 
                                     % subject, 1 to run the preprocessing,  
                                     % 1 to run the 1st level model
                                     
idx=0; % counter variable used to index the jobs structure 
                                     
% initialise job structure:
clear jobs;
jobs(1:njob)=deal(struct('task',[],'n_return_values',[],'input_args',[],'depends_on',[]));


% fill in jobs structure
   % loop through tasks    
for s=1:size(subjects,2) % loop through subjects
    for t=1:size(tasks,2) % loop through tasks for each subject  
        idx=idx+1;
        jobs(idx).task=str2func(tasks{t}); % create a function handle for the current task
        jobs(idx).n_return_values=0;
        jobs(idx).input_args={subjects{s},anadir};
        
        % include information about job dependencies. In this example, the
        % first task (pre-processing) doesn't depend on any other jobs
        % having completed, but the second job (1st level model) can only
        % run after preprocessing has been done:
        if t==1 % preprocessing
        % no dependencies - the depends_on field can be empty or zero
            jobs(idx).depends_on=0;
        elseif t==2 % first level model
        % depends on the previous job id (the pre-processing job for the
        % same subject) having completed
            jobs(idx).depends_on=idx-1;
        end        
    end
end

% NB - in some cases, e.g. if you want to run different tasks for each
% subject, or if you just want to run several unrelated tasks, you might
% have to create the jobs array manually, e.g.

% jobs(1).task = @job_script_1;
% jobs(1).n_return_values = 0;
% jobs(1).input_args = {'foo','bar'};

% jobs(2).task = @job_script_2;
% jobs(2).n_return_values = 0;
% jobs(2).input_args = {'dum'};

% etc



%% ========================================================================
% Create a scheduler object. This allows Matlab to submit jobs to the
% Torque scheduler used on the CBU cluster.

clear scheduler;
scheduler=cbu_scheduler();

% The cbu_scheduler function can be used to configure the scheduler with
% various settings. Called without any input arguments, it will return a
% scheduler configured with a default set of parameters (compute queue, 12
% worker nodes, machines with at least 48Gb RAM). The function will also
% recognise various named configurations, e.g.
%   scheduler = cbu_scheduler('gpu')
%   scheduler = cbu_scheduler('large_compute')
% It will also support custom configurations, e.g.
%   scheduler = cbu_scheduler('custom',{'compute',16,96}))
% For more information, see the cbu_scheduler help section ("help
% cbu_scheduler")


%% ========================================================================
% Submit the jobs to the cbu cluster
cbu_qsub(jobs,scheduler,workerpath);

% (The scheduler argument is optional - if omitted, cbu_qsub will configure
% a scheduler object with a default set of paramters). 
%
% further examples:
%
% cbu_qsub(jobs);
%           - submit jobs to a scheduler configured with default parameters
%
% cbu_qsub(jobs,cbu_scheduler('gpu'))
%           - submit jobs to the gpu queue
%
% cbu_qsub(jobs,cbu_scheduler('custom',{'compute',16,64}))
%           - submit jobs to the compute queue, requesting 16 worker nodes
%           and machines with at least 64Gb of RAM.
