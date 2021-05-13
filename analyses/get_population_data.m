function [population, datalist] = get_population_data(recordinglist, type, path_data, time_window, trial_group)
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
% datalist : table
%     table recording info for each recording/unit
% 


%% specify paths
switch type
    case {'pupil_windows','pupil_timeseries'}
        
        % change folder to load results of analyses
        path_data = regexprep(path_data, 'processed', 'analysed');
        path_data = fullfile(path_data,'pupil_drug_modulation');
        
    case {'spike_rate_summary','spike_rate_ANOVA','spike_rate_ROC','RT_drug_modulation'}
        
        % change folder to load results of analyses
        path_data = regexprep(path_data, 'processed', 'analysed');
        path_data = fullfile(path_data,type);
        
    case 'spike_rate_PSTH'
        path_data = regexprep(path_data, 'processed', 'analysed');
        path_data = fullfile(path_data,'spike_rate_summary');
        
    otherwise
        
end

%% specify filenames
switch type
    case 'pupil_windows'
        filename = sprintf('pupil_window_%s.mat', time_window);
        
    case 'RT_drug_modulation'
        filename = 'RT.mat';
        
    case 'pupil_timeseries'
        filename = sprintf('pupil_timeseries_%s.mat', time_window);

    case 'spike_rate_PSTH'
        filename = sprintf('rate_PSTH.mat');

    case 'spike_rate_summary'
        filename = sprintf('rate_summary_%s_%s.mat', time_window, trial_group);
        
    case 'spike_rate_ANOVA'
        filename = sprintf('rate_ANOVA_%s.mat', time_window);
        
    case 'spike_rate_ROC'
        filename = sprintf('rate_ROC_%s_%s.mat', time_window, trial_group);
   
    case 'waveform'
        filename = 'unit.mat';
end



%% load and concatenate data
datalist = [];
total_unit = 0;
for irec = 1:height(recordinglist)
    
    % define info/paths/filenames for this recording
    recinfo = recordinglist(irec,:);
    tmp_path = fullfile(path_data, recinfo.Subject, recinfo.Date);
    loadfilename = fullfile(tmp_path, filename);
    
    if ~exist(loadfilename,'file')
        warning('file does not exist: %s',loadfilename)
        continue
    end
    loaddata = load(loadfilename);
        
    % collect data
    switch type       
        
        case 'pupil_timeseries'
            
            if irec==1
                population.pupil = [];
                population.time = loaddata.timestamps;
            end
            population.pupil = cat(1, population.pupil, loaddata.pupil_timeseries);

        case 'pupil_windows'
            
            if irec==1
                population.pupil = [];
            end
            population.pupil = cat(1, population.pupil, loaddata.pupil_cond);
                    
        case 'RT_drug_modulation'
            if irec==1
                population.RT = [];
            end
            population.RT = cat(1, population.RT, loaddata.RT);
            
        case 'spike_rate_PSTH'
            num_unit = size(loaddata.PSTH.(time_window).samples,1);

            if irec==1
                population.maxhist = [];
                population.samples = [];
                population.time = loaddata.PSTH.(time_window).time;
            end
                        
            population.maxhist   = cat(1, population.maxhist, loaddata.maxhist);
            population.samples   = cat(1, population.samples, loaddata.PSTH.(time_window).samples);

        case 'spike_rate_summary'
            
            num_unit = size(loaddata.rate_cond,1);
            if irec==1
                population.rate = [];
                population.FF = [];
                population.lm.coef = [];
                population.lm.t = [];
                population.lm.p = [];
                
            end
            
            population.rate = cat(1, population.rate, loaddata.rate_cond);
            population.FF   = cat(1, population.FF, loaddata.FF_cond);
            
            population.lm.coef  = cat(1, population.lm.coef, loaddata.lm.coef);
            population.lm.t     = cat(1, population.lm.t, loaddata.lm.tStat);
            population.lm.p     = cat(1, population.lm.p, loaddata.lm.pValue);

        case 'spike_rate_ANOVA'
            
            num_unit = size(loaddata.p_anova,1);
            if irec==1
                population.selectivity = [];
            end

            p_tmp = [table(loaddata.p_stim, 'VariableNames', {'stim'}), loaddata.p_anova];
            
            population.selectivity = cat(1, population.selectivity, p_tmp);
               
        case 'spike_rate_ROC'
            
            num_unit = size(loaddata.roc_attend,1);

            if irec==1
                population.roc_attend = [];
                population.roc_drug = [];
                population.mi_attend = [];
                population.mi_drug = [];
                population.gain = [];
            end
            
            population.roc_attend   = cat(1, population.roc_attend, loaddata.roc_attend);
            population.roc_drug     = cat(1, population.roc_drug, loaddata.roc_drug);
            population.mi_attend    = cat(1, population.mi_attend, loaddata.mi_attend);
            population.mi_drug      = cat(1, population.mi_drug, loaddata.mi_drug);
            population.gain         = cat(1, population.gain, loaddata.gain);
            
            
        case 'waveform'
            
            num_unit = size(loaddata.waveform.average_spike,1);
            
            if irec==1
                population.waveform = [];
                population.peak_to_trough_time = [];
                
                population.time = loaddata.waveform.time;
            end
            
            population.peak_to_trough_time = cat(1, population.peak_to_trough_time, loaddata.waveform.peak_to_trough_time);
            population.waveform = cat(1, population.waveform, loaddata.waveform.average_spike);
                        
    end
    
    switch type
        case {'pupil_windows','pupil_timeseries','RT_drug_modulation'}
            datalist = [datalist; recinfo];
            
        otherwise
            % concatenate unitlist from recinfo
            datalist = [datalist; [repmat(recinfo, [num_unit,1]), table((1:num_unit)', 'VariableNames', {'unit'}) ]];
    end
    
    if exist('num_unit','var')
        total_unit = total_unit+num_unit;
    end
end


% append table
switch type
    case 'spike_rate_ANOVA'
        
        % replace column names
        varnames = population.selectivity.Properties.VariableNames;
        varnames = cellfun(@(x) (['s_' x]), varnames, 'UniformOutput', false);
        population.selectivity.Properties.VariableNames = varnames;
        
        % add unit column
        population.selectivity.unit = categorical((1:total_unit)', 'Ordinal', false); 

end
