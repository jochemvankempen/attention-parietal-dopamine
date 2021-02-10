function [trialdata, unit] = remove_excluded_trials(trialdata, unit)
% [trialdata, unit] = remove_excluded_trials(trialdata, unit)
%
% remove trials from both trialdata and unit. Trials in the field
% `trialdata.exclude_trial` will be removed.
%

event_fields = fields(unit);
event_fields = event_fields(contains(event_fields, 'Align'));

for ievent = 1:length(event_fields)
    unit.(event_fields{ievent})(:,[trialdata.exclude_trial]) = [];
end
trialdata([trialdata.exclude_trial]) = [];
