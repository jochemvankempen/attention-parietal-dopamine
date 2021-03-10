function [trialdata, unit, idx_exclude_trials] = remove_excluded_trials(trialdata, unit, idx_unit)
% [trialdata, unit] = remove_excluded_trials(trialdata, unit, idx_unit)
%
% remove trials from both trialdata and unit. Trials in the field
% `unit.idx_include_trials` will be removed.
%

event_fields = fields(unit);
event_fields = event_fields(contains(event_fields, 'Align'));

idx_exclude_trials = ...
    ~unit.idx_include_trials(idx_unit,:) ...
    | [trialdata.exclude];

for ievent = 1:length(event_fields)
    unit.(event_fields{ievent})(:,idx_exclude_trials) = [];
end
trialdata(idx_exclude_trials) = [];
