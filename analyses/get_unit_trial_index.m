function trial_index = get_unit_trial_index(unit, unit_index, event)
% trial_index = get_unit_trial_index(unit, unit_index, event)
%
% Select the trial window that has spike for a particular unit
%
% Parameters
% ----------
% unit : struct
%     struct with spiking activity aligned to events listed in fields
%     '...Align'. Fields are of size num_unit x num_trial.
% unit_index : double
%     double indicating which unit to use
% event : string (optional)
%     string indicating which event to use, default is StimAlign.
% 
% Returns
% -------
% trial_index : array of bool
%     booleans indicating trial inclusion
%


if nargin<3
    event = 'StimAlign';
end

[num_unit, num_trial] = size(unit.(event));

spike_count = cellfun(@(x) length(x), unit.(event)(unit_index,:), 'UniformOutput',true);

first_spike = find(spike_count~=0, 1, 'first');
last_spike = find(spike_count~=0, 1, 'last');

trial_index = false(num_trial,1);
trial_index(first_spike:last_spike) = true;
