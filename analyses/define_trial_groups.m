function trialdata = define_trial_groups(trialdata, trial_groups, type)
% trialdata = define_trial_groups(trialdata, trial_groups, type)
%
% Assigns trials to groups based on type
%
% Parameters
% ----------
% trialdata : struct
%     struct with metadata about the trials
% trial_groups : array
%     array of size (num_groups x 2), with columns indicating first and
%     last index of trial groups
% type : string
%     string indicating which trial grouping to apply
% 

% check input, get dimensions
[num_groups, windowsize] = size(trial_groups);
num_trials = length(trialdata);
assert(windowsize==2, 'windowsize!=2, window requires start and end index')


% loop over groups
assigned_group = NaN(num_trials,1);
for igroup = 1:num_groups
    switch type
        case 'block_change'
            block_change = [trialdata.block_change];
            trial_idx = block_change>=trial_groups(igroup,1) & block_change<=trial_groups(igroup,2);
            assigned_group(trial_idx) = igroup;
            
            if igroup==1
                trial_idx = block_change<trial_groups(igroup,1);
                assigned_group(trial_idx) = -1;
            end
    end
    
    if igroup==1
        group_label = sprintf('%d-%d',trial_groups(igroup,1),trial_groups(igroup,2));
    else
        group_label = [group_label '_' sprintf('%d-%d',trial_groups(igroup,1),trial_groups(igroup,2))];
    end
end

% deal
assigned_group = num2cell(assigned_group);
[trialdata.group] = assigned_group{:};
clear assigned_group

% define trialgroup label
[trialdata.group_label] = deal(sprintf('%s(%s)',type, group_label));

% check output
assigned_group = [trialdata.group];
assert(~any(isnan(assigned_group)), 'Not all trials have been assigned a group')