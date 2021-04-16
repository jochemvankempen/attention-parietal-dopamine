function batch_analyses()
% 

%%
% addpath(genpath('C:\Jochem\repositories/attention-parietal-dopamine'))
addpath(genpath('/Users/jochemvankempen/repositories/attention-parietal-dopamine'))
%%

subjects = {'J','W'};% dopamine data
% subjects = {'S'}; % control data

%% paths

% path_data = ['C:\Jochem\Gratc_PPC_DA_new\data\processed'];
path_data = ['/Users/jochemvankempen/NCL/gratc_DA/processed'];

%% addpaths

addpath(genpath(''))

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


