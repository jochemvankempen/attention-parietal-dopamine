function [population, unitlist] = get_population_data(recordinglist, type, path_data, time_window)
% population_data = get_population_data(recordinglist, type, path_data)
%
% Collect and concatenate data from individual recordings into population
% resuls. 
% 
% Parameters
% ----------
% recordinglist : table
%     table with recording info
% type : string
%     string specifying which data to load
% path_data : string
%     string with the location of the datafiles
% time_window : string
%     string specifying which time window to load
%
% Returns
% -------
% population : array/struct/table
%     array/struct/table with concatenated data for the population
% unitlist : table
%     table recording info for each unit
% 


%% specify paths
switch type
    case {'spike_rate_summary','spike_rate_ANOVA','spike_rate_ROC'}
        
        % change folder to load results of analyses
        path_data = regexprep(path_data, 'processed', 'analysed');
        path_data = fullfile(path_data,type);
        
    otherwise
        
end

%% specify filenames
switch type
    case 'spike_rate_summary'
        filename = sprintf('rate_summary_%s.mat', time_window);
        
    case 'spike_rate_ANOVA'
        filename = sprintf('rate_ANOVA_%s.mat', time_window);
        
    case 'spike_rate_ROC'
        filename = sprintf('rate_ROC_%s.mat', time_window);
   
    case 'waveform'
        filename = 'unit.mat';
end



%% load and concatenate data
unitlist = [];
for irec = 1:height(recordinglist)
    
    % define info/paths/filenames for this recording
    recinfo = recordinglist(irec,:);
    tmp_path = fullfile(path_data, recinfo.Subject, recinfo.Date);
    loadfilename = fullfile(tmp_path, filename);
    loaddata = load(loadfilename);
        
    % collect data
    switch type
        case 'spike_rate_summary'
            
            num_unit = size(loaddata.rate_cond,1);
            if irec==1
                population.rate = [];
                population.FF = [];
            end
            
            population.rate = cat(1, population.rate, loaddata.rate_cond);
            population.FF   = cat(1, population.FF, loaddata.FF_cond);
            
            
        case 'spike_rate_ANOVA'
            
            num_unit = size(loaddata.p_anova,1);
            if irec==1
                population.p_anova = [];
            end
            
            population.p_anova = cat(1, population.p_anova, loaddata.p_anova);
               
        case 'spike_rate_ROC'
            
            num_unit = size(loaddata.roc_attend,1);

            if irec==1
                population.roc_attend = [];
                population.roc_drug = [];
                population.mi_attend = [];
                population.mi_drug = [];
            end
            
            population.roc_attend   = cat(1, population.roc_attend, loaddata.roc_attend);
            population.roc_drug     = cat(1, population.roc_drug, loaddata.roc_drug);
            population.mi_attend    = cat(1, population.mi_attend, loaddata.mi_attend);
            population.mi_drug      = cat(1, population.mi_drug, loaddata.mi_drug);
            
            
        case 'waveform'
            
            num_unit = size(loaddata.waveform.average_spike,1);
            
            if irec==1
                population.waveform = [];
                population.peak_to_trough_time = [];
                
                population.waveform_time = loaddata.waveform.time;
            end
            
            population.peak_to_trough_time = cat(1, population.peak_to_trough_time, loaddata.waveform.peak_to_trough_time);
            population.waveform = cat(1, population.waveform, loaddata.waveform.average_spike);
                        
    end
    
    % concatenate unitlist from recinfo
    unitlist = [unitlist; repmat(recinfo, [num_unit,1])];

end







