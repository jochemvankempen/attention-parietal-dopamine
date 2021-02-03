function spike_rate_summary(recinfo, trialdata, unit, time_windows, path_target)
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
% path_target : string
%     path where to store results
% 

% check input 
assert(height(recinfo)==1, 'height recInfo ~= 1')

% get dimensions
[num_unit, num_trials] = size(unit.StimAlign);
conditions = unique([trialdata.cond_num]);

% compute firing rate
rate = spike_rate(unit, time_windows);

% get alignments
timewin_fields = fields(rate);

% define path to store result
path_target = fullfile(path_target, 'spike_rate_summary', recinfo.Subject, recinfo.Date);
if ~isfolder(path_target)
    mkdir(path_target)
end

% loop over time windows
for itw = 1:length(timewin_fields)
    
    [rate_cond, FF_cond] = deal(NaN(num_unit, 2, length(conditions)));% unit, drug, cond
    
    for iunit = 1:num_unit
        
        % select only trial window for which this unit has spikes
        trial_index = get_unit_trial_index(unit, iunit);
        tmp_rate = rate.(timewin_fields{itw})(iunit,trial_index)';
        tmp_trialdata =  trialdata(trial_index);
        
        idx_cond = [tmp_trialdata.cond_num]';
        idx_drug = [tmp_trialdata.drug]' + 1;

        % rate across condition
        for idrug = 1:2
            for icond = 1:length(conditions)
                
                % get trial indices
                trial_index = (...
                    idx_drug==idrug) ...
                    & idx_cond==icond;
                
                % rate
                [rate_cond(iunit,idrug,icond)] = mean(tmp_rate(trial_index));
                
                % FF
                [FF_cond(iunit,idrug,icond)] = var(tmp_rate(trial_index)) / mean(tmp_rate(trial_index));
        
            end
        end
        clear tmp*
    end
        
    savefilename = fullfile(path_target, sprintf('rate_ROC_%s.mat', timewin_fields{itw}));
    save(savefilename, 'rate_*', 'FF_*', 'time_windows');
    
end