function spike_rate_ROC(recinfo, trialdata, unit, time_windows, path_target)
% spike_rate_ROC(recinfo, spike_rate, path_target)
%
% Computes attention/drug ROC over each field in time_windows
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
num_group = length(unique([trialdata.group]));

% compute firing rate
rate = spike_rate(unit, time_windows, 'rate');

% get alignments
timewin_fields = fields(rate);

% define path to store result
path_target = fullfile(path_target, 'spike_rate_ROC', recinfo.Subject, recinfo.Date);
if ~isfolder(path_target)
    mkdir(path_target)
end

% set task specific param
switch recinfo.Task
    case 'gratc'
        cond_compare = [...
            1 2;
            1 3;
            3 2;
            4 5;
            4 6;
            6 5;];
               
    case 'msacc'
        error('not implemented yet')
end

% loop over time windows
for itw = 1:length(timewin_fields)
    
    [roc_attend, mi_attend] = deal(NaN(num_unit, 2, size(cond_compare,1)), num_group); % unit, drug, cond, group
    [roc_drug, mi_drug] = deal(NaN(num_unit, size(cond_compare,1)), num_group); % unit, cond, group

    for iunit = 1:num_unit
        
        % select only trial window for which this unit has spikes
        [tmp_trialdata, tmp_unit, idx_exclude_trials] = remove_excluded_trials(trialdata, unit, iunit);
        tmp_rate = rate.(timewin_fields{itw})(iunit,~idx_exclude_trials)';
        
        idx_cond = [tmp_trialdata.cond_num]';
        idx_drug = [tmp_trialdata.drug]' + 1;
        idx_group = [tmp_trialdata.group]';
        
        % attention ROC
        for idrug = 1:2
            for icond = 1:length(cond_compare)
                for igroup = 1:num_group
                    
                    % get trial indices
                    trial_index = ...
                        idx_drug==idrug ...
                        & (idx_cond==cond_compare(icond,1) | idx_cond==cond_compare(icond,2)) ...
                        & idx_group==igroup;
                    
                    % define trial classes
                    class = NaN(length(trial_index),1);
                    class(idx_cond==cond_compare(icond,1)) = 1;
                    class(idx_cond==cond_compare(icond,2)) = 0;
                    
                    % ROC analysis
                    [roc_attend(iunit,idrug,icond,igroup),~,~] = ROC_area(tmp_rate(trial_index), class(trial_index));
                    
                    % attMI = (attRF - attAway) / (attRF + attAway)
                    tmp = ...
                        ( mean(tmp_rate(trial_index & class==1)) - mean(tmp_rate(trial_index & class==0)) ) / ...
                        ( mean(tmp_rate(trial_index & class==1)) + mean(tmp_rate(trial_index & class==0)) );
                    mi_attend(iunit,idrug,icond,igroup) = tmp;
                    
                end
            end
        end
        
        % drug ROC
        for icond = 1:length(unique(idx_cond))
            for igroup = 1:num_group
                
                % get trial indices
                trial_index = ...
                    idx_cond==icond ...
                    & idx_group==igroup;
                
                % define trial classes
                class = NaN(length(trial_index),1);
                class(idx_drug==1) = 0;
                class(idx_drug==2) = 1;
                
                % ROC analysis
                [roc_drug(iunit,icond,igroup),~,~] = ROC_area(tmp_rate(trial_index), class(trial_index));
                
                % drugMI = (drug - no drug) / (drug + no drug)
                tmp = ...
                    ( mean(tmp_rate(trial_index & class==1)) - mean(tmp_rate(trial_index & class==0)) ) / ...
                    ( mean(tmp_rate(trial_index & class==1)) + mean(tmp_rate(trial_index & class==0)) );
                
                mi_drug(iunit,icond,igroup) = tmp;
            end
        end
        
        clear tmp*
    end
        
    savefilename = fullfile(path_target, sprintf('rate_ROC_%s_%s.mat', timewin_fields{itw}, trialdata(1).group_label));
    save(savefilename, 'roc_*', 'mi_*', 'time_windows');
    
end