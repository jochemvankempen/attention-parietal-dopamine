function perform_analyses(recinfo, Job, path_data)
% perform_analyses(recinfo, Job, path_data)
%
% Wrapper function from which analyses/Jobs are called
%
% Parameters
% ----------
% recinfo : table
%     table with single row (one recording)
% Job : struct
%     structure with fields refering to individual analyses. Each field
%     is a boolean indicating whether the analysis will be executed
% path_data : string
%     string with the location of the datafiles
%

%% check input

assert(height(recinfo)==1, 'height recInfo ~= 1')

%% param

% define analysis timewindows for statistical analyses
time_windows.baseline = {'Stim', [-200 0]};
time_windows.stim = {'Stim', [0 400]};
time_windows.cue = {'Cue', [0 400]};
time_windows.pre_cue = {'Cue', [-200 0]};
time_windows.dim = {'Dim1', [-500 0]};

%% load data

loadfilename = fullfile(path_data, recinfo.Subject, recinfo.Date, 'unit.mat');
unit = load(loadfilename);
loadfilename = fullfile(path_data, recinfo.Subject, recinfo.Date, 'trialdata.mat');
load(loadfilename);

path_target = regexprep(path_data, 'processed', 'analysed'); % path where results are stored

%% remove trials 

event_fields = fields(unit);
event_fields = event_fields(contains(event_fields, 'Align'));

for ievent = 1:length(event_fields)
    unit.(event_fields{ievent})(:,[trialdata.exclude_trial]) = [];
end
trialdata([trialdata.exclude_trial]) = [];


%% mean rate & FF

if Job.spike_rate_summary
    spike_rate_summary(recinfo, trialdata, unit, time_windows, path_target)
end

%% perform anova across spike rates 

if Job.spike_rate_ANOVA
    spike_rate_ANOVA(recinfo, trialdata, unit, time_windows, path_target)
end

%% attention/drug ROC 

if Job.spike_rate_ROC
    spike_rate_ROC(recinfo, trialdata, unit, time_windows, path_target)
end
















