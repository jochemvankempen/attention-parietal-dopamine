function spike_rate = spike_rate(unit, time_windows)
% spike_rate = spike_rate(unit, time_windows)
%
% Computes the firing rate within time windows defined in  `time_windows`
%
% Parameters
% ----------
% unit : struct
%     
% time_windows : struct
%     
%
% Returns
% -------
% spike_rate : struct
%     struct with firing rate per time window
%


time_fields = fields(time_windows);

for itime = 1:length(time_fields)
    
    event = time_windows.(time_fields{itime}){1};
    time_window = time_windows.(time_fields{itime}){2};
    
    tmp_count = cellfun(@(x) length(find(x>time_window(1) & x<time_window(2))), unit.([event 'Align']), 'UniformOutput',true);
    tmp_rate = tmp_count / (diff(time_window)/1000);
    
    spike_rate.(time_fields{itime}) = tmp_rate;
end
