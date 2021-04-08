function pupil_drug_modulation(recinfo, trialdata, time_windows, timeseries_windows, path_target)
% pupil_drug_modulation(recinfo, trialdata, time_windows, timeseries_windows, path_target)
%
% Computes average pupil diameter across conditions
%
% Parameters
% ----------
% recinfo : table
%     table with single row (one recording)
% trialdata : struct
%     struct with metadata about the trials
% time_windows : struct
%     struct with time windows  to evaluate
% PSTH_windows : struct
%     struct with time windows for which to extract PSTH
% path_target : string
%     path where to store results
% 

% check input 
assert(height(recinfo)==1, 'height recInfo ~= 1')

% check whether pupil data exists
trialdata_fields = fields(trialdata);
if ~any(strcmpi(trialdata_fields, 'pupil'))
    return
end

% remove excluded trials
trialdata([trialdata.exclude]) = [];

% get metadata
num_trial = length(trialdata);
conditions = unique([trialdata.cond_num]);
idx_cond = [trialdata.cond_num]';
idx_drug = [trialdata.drug]' + 1;


% get alignments
timewin_fields = fields(time_windows);
timeseries_fields = fields(timeseries_windows);

% define path to store result
path_target = fullfile(path_target, 'pupil_drug_modulation', recinfo.Subject, recinfo.Date);
if ~isfolder(path_target)
    mkdir(path_target)
end

% find normalisation factor
% samples_min = max(arrayfun(@(x) min(x.pupil.samples), trialdata, 'UniformOutput', true));
% samples_max = max(arrayfun(@(x) max(x.pupil.samples), trialdata, 'UniformOutput', true));
samples_mean = mean(arrayfun(@(x) nanmean(x.pupil.samples), trialdata, 'UniformOutput', true));
samples_std = mean(arrayfun(@(x) nanstd(x.pupil.samples), trialdata, 'UniformOutput', true));

% extract average pupil diameter in given time windows
% ----------------------------------------------------
for itrial = 1:num_trial

    time = trialdata(itrial).pupil.timestamps;
%     samples = trialdata(itrial).pupil.samples / samples_max;
    samples = (trialdata(itrial).pupil.samples - samples_mean)/samples_std; % z-score
    
    for itw = 1:length(timewin_fields)
        
        % which event to align to
        event = time_windows.(timewin_fields{itw}){1};
        event_time = trialdata(itrial).pupil.events.(event);
        
        % get timestamps
        time_event = time - event_time; % align to event
        
        % which window to average in
        window = time_windows.(timewin_fields{itw}){2};

        % find time indices for window
        idx_time = dsearchn(time_event, window(:));
        
        % average pupil in time window
        pupil_window.(timewin_fields{itw})(itrial,1) = nanmean(samples(idx_time(1):idx_time(2)));
    end
end

% baseline subtraction, go over time windows again to subtract baseline
for itw = 1:length(timewin_fields)
    
    if strcmpi(timewin_fields{itw}, 'baseline')
        continue
    end
    
    pupil_window.(timewin_fields{itw}) = pupil_window.(timewin_fields{itw}) - pupil_window.baseline;
end

% store average for each condition
for itw = 1:length(timewin_fields)
    [pupil_cond] = deal(NaN(1,2,length(conditions)));% drug, cond
    
    for idrug = 1:2
        for icond = 1:length(conditions)
            
            % get trial indices
            trial_index = (...
                idx_drug==idrug) ...
                & idx_cond==icond;
            
            % pupil
            [pupil_cond(1,idrug,icond)] = mean(pupil_window.(timewin_fields{itw})(trial_index));
        end
    end
    
    savefilename = fullfile(path_target, sprintf('pupil_window_%s.mat', timewin_fields{itw}));
    save(savefilename, 'pupil_cond', 'time_windows');

end

% extract aligned time series
% ---------------------------
for itrial = 1:num_trial
    
    time = trialdata(itrial).pupil.timestamps;
    % samples = trialdata(itrial).pupil.samples / samples_max;    
    samples = (trialdata(itrial).pupil.samples - samples_mean)/samples_std; % z-score
    samples = samples - pupil_window.baseline(itrial); % subtract baseline
    
        
    for itw = 1:length(timeseries_fields)
        
        % which event to align to
        event = timeseries_windows.(timeseries_fields{itw}){1};
        event_time = trialdata(itrial).pupil.events.(event);
        
        % get timestamps
        time_event = time - event_time; % align to event
        
        % which window to average in
        alignment = timeseries_windows.(timeseries_fields{itw}){2};

        % find time indices for window
        idx_time = dsearchn(time_event, alignment(:));
        
        % number of samples to extract
        num_samples = timeseries_windows.(timeseries_fields{itw}){3};
        
        if itrial==1        
            % init
            pupil_ts.(timeseries_fields{itw}) = NaN(num_trial, num_samples);
            
            pupil_time.([timeseries_fields{itw}]) = time_event(idx_time(1):idx_time(1)+num_samples-1)';
        end
        
        % average pupil in time window
        pupil_ts.(timeseries_fields{itw})(itrial,:) = samples(idx_time(1):(idx_time(1)+num_samples-1));
        
    end
end


% store average for each condition
for itw = 1:length(timeseries_fields)
    timestamps = pupil_time.([timeseries_fields{itw}]);
    
    [pupil_timeseries] = deal(NaN(1,2,length(conditions),size(pupil_ts.(timeseries_fields{itw}),2)));% drug, cond, time

    for idrug = 1:2
        for icond = 1:length(conditions)
            
            % get trial indices
            trial_index = (...
                idx_drug==idrug) ...
                & idx_cond==icond;
            
            % pupil
            pupil_timeseries(1,idrug,icond,:) = nanmean(pupil_ts.(timeseries_fields{itw})(trial_index,:),1);
        end
    end

    savefilename = fullfile(path_target, sprintf('pupil_timeseries_%s.mat', timeseries_fields{itw}));
    save(savefilename, 'timestamps', 'pupil_timeseries', 'timeseries_windows');

end


