function batch_population()


%% param
subjects = {'J','W'};
path_data = ['/Users/jochemvankempen/NCL/gratc_DA/processed'];

%% get recordinglist
recordinglist = get_recordinglist(subjects, path_data);

%% 

% rate_summary = get_population_data(recordinglist, 'spike_rate_summary', path_data, 'dim');
% rate_ANOVA = get_population_data(recordinglist, 'spike_rate_ANOVA', path_data, 'dim');
% rate_ROC = get_population_data(recordinglist, 'spike_rate_ROC', path_data, 'dim');
waveform = get_population_data(recordinglist, 'waveform', path_data);

% selectivity = get_unit_selectivity(rate_ANOVA.p_anova, 'att+dru')

keyboard