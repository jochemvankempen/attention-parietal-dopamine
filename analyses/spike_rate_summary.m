function spike_rate_summary(recinfo, trialdata, unit, time_windows, PSTH_windows, path_target)
% spike_rate_summary(recinfo, spike_rate, path_target)
%
% Computes average rate and fano factor for each condition over each field in time_windows
%
% Parameters
% ----------
% recinfo : table
%     table with single row (one recording)
% trialdata : struct
%     struct with metadata about the trials
% unit : struct
%     struct with spike times aligned to events
% time_windows : struct
%     struct with time windows  to evaluate
% PSTH_windows : struct
%     struct with time windows for which to extract PSTH
% path_target : string
%     path where to store results
% 

% check input 
assert(height(recinfo)==1, 'height recInfo ~= 1')

% get dimensions
[num_unit, num_trials] = size(unit.StimAlign);
conditions = unique([trialdata.cond_num]);

% compute firing rate
rate = spike_rate(unit, time_windows, 'rate');

% get alignments
timewin_fields = fields(rate);
PSTH_fields = fields(PSTH_windows);

% define path to store result
path_target = fullfile(path_target, 'spike_rate_summary', recinfo.Subject, recinfo.Date);
if ~isfolder(path_target)
    mkdir(path_target)
end

% compute histogram
binsize = 30;
maxhist = zeros(num_unit,1); % normalisation factor
for itw = 1:length(PSTH_fields)
        
    PSTH_event = PSTH_windows.(PSTH_fields{itw}){1};
    time_window = PSTH_windows.(PSTH_fields{itw}){2};
    
    for iunit = 1:num_unit
        
        % select only trial window for which this unit has spikes
        [tmp_trialdata, tmp_unit] = remove_excluded_trials(trialdata, unit, iunit);

        % compute histogram
        [hist_spike, hist_time] = spike_convolute_gaussian(tmp_unit.([PSTH_event 'Align']), time_window, iunit, binsize);
        
        % init store
        if iunit==1
            PSTH.(PSTH_fields{itw}).samples = NaN(num_unit, 2, length(conditions), length(hist_time)); % unit, drug, conditions, time
            PSTH.(PSTH_fields{itw}).time = hist_time;
        end
        
        % rate across condition
        idx_cond = [tmp_trialdata.cond_num]';
        idx_drug = [tmp_trialdata.drug]' + 1;
        
        for idrug = 1:2
            for icond = 1:length(conditions)
                
                % get trial indices
                idx_trial = (...
                    idx_drug==idrug) ...
                    & idx_cond==icond;
                
                % store
                PSTH.(PSTH_fields{itw}).samples(iunit,idrug,icond,:) = mean(hist_spike(idx_trial,:),1);
                
                if maxhist(iunit)<max(PSTH.(PSTH_fields{itw}).samples(iunit,idrug,icond,:))
                    maxhist(iunit) = max(PSTH.(PSTH_fields{itw}).samples(iunit,idrug,icond,:));
                end
            end
        end
    end
end

% normalisation
for itw = 1:length(PSTH_fields)
    for iunit = 1:num_unit
        PSTH.(PSTH_fields{itw}).samples(iunit,:,:,:) = PSTH.(PSTH_fields{itw}).samples(iunit,:,:,:) / maxhist(iunit) * 100; 
    end
end
savefilename = fullfile(path_target, sprintf('rate_PSTH.mat'));
save(savefilename, 'PSTH', 'PSTH_windows', 'maxhist');


% loop over time windows
for itw = 1:length(timewin_fields)   
    
    % init
    [rate_cond, FF_cond] = deal(NaN(num_unit, 2, length(conditions)));% unit, drug, cond
    
    for iunit = 1:num_unit
                
        % select only trial window for which this unit has spikes
        [tmp_trialdata, tmp_unit, idx_exclude_trials] = remove_excluded_trials(trialdata, unit, iunit);
        tmp_rate = rate.(timewin_fields{itw})(iunit,~idx_exclude_trials)';
        
        idx_cond = [tmp_trialdata.cond_num]';
        idx_drug = [tmp_trialdata.drug]' + 1;

        % rate across condition
        for idrug = 1:2
            for icond = 1:length(conditions)
                
                % get trial indices
                idx_trial = (...
                    idx_drug==idrug) ...
                    & idx_cond==icond;
                
                % rate
                [rate_cond(iunit,idrug,icond)] = mean(tmp_rate(idx_trial));
                
                % FF
                [FF_cond(iunit,idrug,icond)] = var(tmp_rate(idx_trial)) / mean(tmp_rate(idx_trial));
        
            end
        end
        clear tmp*
    end
        
    savefilename = fullfile(path_target, sprintf('rate_summary_%s.mat', timewin_fields{itw}));
    save(savefilename, 'rate_*', 'FF_*', 'time_windows');
    
end