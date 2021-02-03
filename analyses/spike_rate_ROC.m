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

% compute firing rate
rate = spike_rate(unit, time_windows);

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
            2 3;
            4 5;
            4 6;
            5 6;];
               
    case 'msacc'
        error('not implemented yet')
end

% loop over time windows
for itw = 1:length(timewin_fields)
    
    [roc_attend, mi_attend] = deal(NaN(num_unit, 2, size(cond_compare,1))); % unit, drug, cond
    [roc_drug, mi_drug] = deal(NaN(num_unit, size(cond_compare,1))); % unit, cond

    for iunit = 1:num_unit
        
        % select only trial window for which this unit has spikes
        trial_index = get_unit_trial_index(unit, iunit);
        tmp_rate = rate.(timewin_fields{itw})(iunit,trial_index)';
        tmp_trialdata =  trialdata(trial_index);
        
        idx_cond = [tmp_trialdata.cond_num]';
        idx_drug = [tmp_trialdata.drug]' + 1;

        % attention ROC
        for idrug = 1:2
            for icond = 1:length(cond_compare)
                
                % get trial indices
                trial_index = (...
                    idx_drug==idrug) ...
                    & (idx_cond==cond_compare(icond,1) | idx_cond==cond_compare(icond,2));
                
                % define trial classes
                class = NaN(length(trial_index),1);
                class(idx_cond==cond_compare(icond,1)) = 1;
                class(idx_cond==cond_compare(icond,2)) = 0;
                
                % ROC analysis
                [roc_attend(iunit,idrug,icond),~,~] = ROC_area(tmp_rate(trial_index), class(trial_index));
        
                % attMI = (attRF - attAway) / (attRF + attAway)
                tmp = ...
                    ( mean(tmp_rate(trial_index & class==1)) - mean(tmp_rate(trial_index & class==0)) ) / ...
                    ( mean(tmp_rate(trial_index & class==1)) + mean(tmp_rate(trial_index & class==0)) );
                mi_attend(iunit,idrug,icond) = tmp;
                
            end
        end
        
        
        % drug ROC
        for icond = 1:length(cond_compare)
            
            % get trial indices
            trial_index = (idx_cond==icond);
            
            % define trial classes
            class = NaN(length(trial_index),1);
            class(idx_drug==1) = 0;
            class(idx_drug==2) = 1;

            % ROC analysis
            [roc_drug(iunit,icond),~,~] = ROC_area(tmp_rate(trial_index), class(trial_index));

            % drugMI = (drug - no drug) / (drug + no drug)
            tmp = ...
                ( mean(tmp_rate(trial_index & class==1)) - mean(tmp_rate(trial_index & class==0)) ) / ...
                ( mean(tmp_rate(trial_index & class==1)) + mean(tmp_rate(trial_index & class==0)) );
            
            mi_drug(iunit,icond) = tmp;
            
        end
        
        clear tmp*
    end
        
    savefilename = fullfile(path_target, sprintf('rate_ROC_%s.mat', timewin_fields{itw}));
    save(savefilename, 'roc_*', 'mi_*', 'time_windows');
    
end