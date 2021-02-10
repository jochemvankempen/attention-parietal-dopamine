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

% define timewindows for further analyses. All the below analyses will be
% run across each of these time windows.
time_windows.baseline = {'Stim', [-200 0]};
time_windows.stim = {'Stim', [0 400]};
time_windows.cue = {'Cue', [0 400]};
time_windows.pre_cue = {'Cue', [-200 0]};
time_windows.dim = {'Dim1', [-500 0]};

PSTH_windows.stim = {'Stim', [-250 1000]};
PSTH_windows.cue = {'Cue', [-250 1000]};
PSTH_windows.dim = {'Dim1', [-1000 250]};

pupil_windows.baseline = {'STIM_ON', [-300 -50]/1000};
pupil_windows.stim = {'STIM_ON', (500 + [-125 125])/1000};
pupil_windows.cue = {'CUE_ON', (500 + [-125 125])/1000};
pupil_windows.dim = {'DIMMING1', [-300 -50]/1000};

pupil_timeseries.stim = {'STIM_ON', -0.5, 501}; % event, timestamp, num_timestamps
pupil_timeseries.cue = {'CUE_ON', -0.5, 501};
pupil_timeseries.dim = {'DIMMING1', -1.0, 301};

%% load data

loadfilename = fullfile(path_data, recinfo.Subject, recinfo.Date, 'unit.mat');
unit = load(loadfilename);
loadfilename = fullfile(path_data, recinfo.Subject, recinfo.Date, 'trialdata.mat');
load(loadfilename);

path_target = regexprep(path_data, 'processed', 'analysed'); % path where results are stored

%% remove trials 

[trialdata, unit] = remove_excluded_trials(trialdata, unit);

%% mean rate & FF

if Job.spike_rate_summary
    spike_rate_summary(recinfo, trialdata, unit, time_windows, PSTH_windows, path_target)
end

%% perform anova across spike rates 

if Job.spike_rate_ANOVA
    spike_rate_ANOVA(recinfo, trialdata, unit, time_windows, path_target)
end

%% attention/drug ROC 

if Job.spike_rate_ROC
    spike_rate_ROC(recinfo, trialdata, unit, time_windows, path_target)
end

%% pupil

if Job.pupil_drug_modulation
    pupil_drug_modulation(recinfo, trialdata, pupil_windows, pupil_timeseries, path_target)
end

%% RT

if Job.RT_drug_modulation
    RT_drug_modulation(recinfo, trialdata, path_target)
end

