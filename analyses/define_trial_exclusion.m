function trialdata = define_trial_exclusion(trialdata, exclusion_criteria, type)
% trialdata = define_trial_exclusion(trialdata, trial_groups, type)
%
% Exclude trials based on type
%
% Parameters
% ----------
% trialdata : struct
%     struct with metadata about the trials
% trial_groups : array
%     array with indices indicating trial in/exclusion thresholds
% type : string
%     string indicating which trial exclusion to apply
% 

% check input, get dimensions
num_trials = length(trialdata);

% define excluded trials
idx_exclude = false(num_trials,1);
switch type
    case 'block_change'

        block_change = [trialdata.block_change]';
        idx_block_change = block_change <= exclusion_criteria;
        
        % exclude first block, as no drug has been applied yet
        idx_block = true(num_trials,1);
        idx_block_exclude = find(diff([trialdata.block_num]),1,'first');
        idx_block(1:idx_block_exclude) = false;

        idx_exclude(idx_block_change & idx_block) = true;
        
end
    

% deal
idx_exclude = num2cell(idx_exclude);
[trialdata.exclude] = idx_exclude{:};
