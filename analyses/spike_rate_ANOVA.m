function spike_rate_ANOVA(recinfo, trialdata, unit, time_windows, path_target)
% spike_rate_ANOVA(recinfo, spike_rate, path_target)
%
% Computes an ANOVA over each field in time_windows
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
count = spike_rate(unit, time_windows, 'count');
rate = spike_rate(unit, time_windows, 'rate');

% get alignments
timewin_fields = fields(count);

% define path to store result
path_target = fullfile(path_target, 'spike_rate_ANOVA', recinfo.Subject, recinfo.Date);
if ~isfolder(path_target)
    mkdir(path_target)
end

% loop over time windows
p_stim = NaN(num_unit,1);
for itw = 1:length(timewin_fields)
    
    % define parameters for ANOVA
    switch recinfo.Task
        case 'gratc'
            attention_text = {'attend RF','attend away'};
            direction_text = {'dir 1','dir 2'};
            drug_text = {'no drug','drug'};
            
            attention = [1 2 2 1 2 2];
            direction = [1 1 1 2 2 2];
            
            group_names{1,1} = 'att';
            group_names{2,1} = 'dru';
            group_names{3,1} = 'dir';
            
            P = zeros(num_unit, 7);
            
%             P = zeros(num_unit, 3);

        case 'msacc'
            condition_text = {'1','2','3','4'};
            drug_text = {'no drug','drug'};
            
            group_names{1,1} = 'cond';
            group_names{2,1} = 'dru';
            
            P = zeros(num_unit, 3);
    end
    
    statstable = cell(num_unit,1);
    stats = cell(num_unit,1);
    terms = cell(num_unit,1);
    for iunit = 1:num_unit
        
        % select only trial window for which this unit has spikes
        [tmp_trialdata, tmp_unit, idx_exclude_trials] = remove_excluded_trials(trialdata, unit, iunit);
        tmp_count = count.(timewin_fields{itw})(iunit,~idx_exclude_trials)';
        
        % define groups for anova
        switch recinfo.Task
            case 'gratc'
                factor.attention = attention([tmp_trialdata.cond_num])';
                factor.drug = [tmp_trialdata.drug]' + 1;
                factor.direction = direction([tmp_trialdata.cond_num])';
                groups = {attention_text(factor.attention)', drug_text(factor.drug)', direction_text(factor.direction)'};
                
%                 groups = {attention_text(factor.attention)', drug_text(factor.drug)'};
                
            case 'msacc'
                factor.condition = [tmp_trialdata.cond_num]';
                factor.drug = [tmp_trialdata.drug]' + 1;
                
                groups = {condition_text(factor.condition)', drug_text(factor.drug)'};
        end
        
        % ANOVA
        [P(iunit,:),statstable{iunit},stats{iunit},terms{iunit}]  = anovan(tmp_count,groups,'full',3,group_names,'off');
        
        % check for significant stimulus response
        switch timewin_fields{itw}
            case 'stim'
                
                idx_baseline = strcmpi(timewin_fields, 'baseline');
                
                rate_baseline = rate.(timewin_fields{idx_baseline})(iunit,~idx_exclude_trials)';
                rate_stim = rate.(timewin_fields{itw})(iunit,~idx_exclude_trials)';
                
                [p_stim(iunit), H, STATS] = signrank(rate_stim-rate_baseline);
                      
%                 3.3030
%     4.7497

        end
        
    end
    p_anova = array2table(P, 'VariableNames', statstable{iunit}(2:(end-2),1));
        
    savefilename = fullfile(path_target, sprintf('rate_ANOVA_%s.mat', timewin_fields{itw}));
    save(savefilename, 'p_*', 'statstable', 'stats', 'terms', 'time_windows');
    
end

% % 
% P =
% 
%     0.0024
%     0.5727
%     0.8328
%     0.5523
%     0.5252
%     0.3083
%     0.6414