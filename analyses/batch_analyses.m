function batch_analyses()
% 

%%
addpath(genpath('M:\Jochem\papers\2021-DA-LIP\attention-parietal-dopamine'))
addpath(genpath('M:\Jochem\papers\2021-DA-LIP\repositories\gain-variability'))
% addpath(genpath('/Users/jochemvankempen/repositories/attention-parietal-dopamine'))
% addpath(genpath('/Users/jochemvankempen/repositories/gain-variability'))
%%

subjects = {'J','W'};% dopamine data
% subjects = {'S'}; % control data

%% paths

path_data = ['M:\Jochem\papers\2021-DA-LIP\data\processed'];
% path_data = ['/Users/jochemvankempen/NCL/gratc_DA/processed'];

%% addpaths

addpath(genpath(''))

%% get list of recordings

recordinglist = get_recordinglist(subjects, path_data);

%% define jobs

Job.spike_rate_summary      = 0;
Job.spike_rate_ANOVA        = 0;
Job.spike_rate_ROC          = 1;
Job.pupil_drug_modulation   = 0;
Job.RT_drug_modulation      = 0;

%% Run jobs

for irec = 1:height(recordinglist)
    
    recinfo = recordinglist(irec,:);
    
    fprintf('Running analyses for: Subject %s, recording %s\n', recinfo.Subject, recinfo.Date)
    fprintf('--------------------------------------------------------\n')
    
    perform_analyses(recinfo, Job, path_data)
        
end


