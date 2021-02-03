function perform_analyses(recinfo, Job, path_data)


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


%% perform anova across spike rates 

if Job.spike_rate_ANOVA
    spike_rate_ANOVA(recinfo, trialdata, unit, time_windows, path_target)
end

%% test attentional ROC 

if Job.spike_rate_ROC
    spike_rate_ROC(recinfo, trialdata, unit, time_windows, path_target)
end

%% Attentional/drug modulation - ROC

% define which conditions to compare
cond_compare = [...
    1 2;
    1 3;
    2 3;
    4 5;
    4 6;
    5 6;];


conditions = unique([trialdata.cond_num]);

nchoosek(conditions, 2)

ROC_area


keyboard


















