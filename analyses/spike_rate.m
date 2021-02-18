function spike_rate = spike_rate(unit, time_windows, type)
% spike_rate = spike_rate(unit, time_windows, type)
%
% Computes the firing count/rate within time windows defined in  `time_windows`
%
% Parameters
% ----------
% unit : struct
%     
% time_windows : struct
%     
% type : string
%     string indicating whether to return spike count or rate (default)
%
% Returns
% -------
% spike_rate : struct
%     struct with firing rate per time window
%

if nargin<3
    type='rate';
end

time_fields = fields(time_windows);

for itime = 1:length(time_fields)
    
    event = time_windows.(time_fields{itime}){1};
    time_window = time_windows.(time_fields{itime}){2};
    
    tmp_count = cellfun(@(x) length(find(x>time_window(1) & x<time_window(2))), unit.([event 'Align']), 'UniformOutput',true);
    switch type
        case 'rate'
            tmp_count = tmp_count * (1000/diff(time_window));
        case 'count'
        otherwise
            error('unknown spike count conversion type')
    end
    spike_rate.(time_fields{itime}) = tmp_count;
end
