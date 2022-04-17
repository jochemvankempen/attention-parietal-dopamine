function batch_analyses()
% run analysis on individual recordings, store output for later population
% level analyses

%% setup repositories, addpath

path_repo = fileparts(mfilename('fullpath'));

% https://github.com/jochemvankempen/attention-parietal dopamine
addpath(genpath(fullfile(path_repo,'/../../attention-parietal-dopamine')))

% https://github.com/jochemvankempen/gain-variability
addpath(genpath(fullfile(path_repo,'/../../gain-variability')))
%%

subjects = {'J','W'};% dopamine data

%% paths

path_data = fullfile(path_repo,'/../../../data/processed');

%% get list of recordings

recordinglist = get_recordinglist(subjects, path_data);

%% define jobs

Job.spike_rate_summary      = 1;
Job.spike_rate_ANOVA        = 1;
Job.spike_rate_ROC          = 1;
Job.pupil_drug_modulation   = 1;
Job.RT_drug_modulation      = 1;

%% Run jobs

for irec = 1:height(recordinglist)
    
    recinfo = recordinglist(irec,:);
    
    fprintf('Running analyses for: Subject %s, recording %s\n', recinfo.Subject, recinfo.Date)
    fprintf('--------------------------------------------------------\n')
    
    perform_analyses(recinfo, Job, path_data)
        
end


