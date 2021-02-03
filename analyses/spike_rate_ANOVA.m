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
rate = spike_rate(unit, time_windows);

% get alignments
timewin_fields = fields(rate);

% define path to store result
path_target = fullfile(path_target, 'spike_rate_ANOVA', recinfo.Subject, recinfo.Date);
if ~isfolder(path_target)
    mkdir(path_target)
end

% loop over time windows
for itw = 1:length(timewin_fields)
    
    % define parameters for ANOVA
    switch recinfo.Task
        case 'gratc'
            attention_text = {'attend RF','attend away'};
            direction_text = {'dir 1','dir 2'};
            drug_text = {'no drug','drug'};
            
            attention = [1 2 2 1 2 2];
            direction = [1 2 1 2 1 2];
            
            group_names{1,1} = 'att';
            group_names{2,1} = 'dru';
            group_names{3,1} = 'dir';
            
            P = zeros(num_unit, 7);
            
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
        trial_index = get_unit_trial_index(unit, iunit);
        tmp_rate = rate.(timewin_fields{itw})(iunit,trial_index)';
        tmp_trialdata =  trialdata(trial_index);
        
        % define groups for anova
        switch recinfo.Task
            case 'gratc'
                factor.attention = attention([tmp_trialdata.cond_num])';
                factor.drug = [tmp_trialdata.drug]' + 1;
                factor.direction = direction([tmp_trialdata.cond_num])';
                
                groups = {attention_text(factor.attention)', drug_text(factor.drug)', direction_text(factor.direction)'};
                
            case 'msacc'
                factor.condition = [tmp_trialdata.cond_num]';
                factor.drug = [tmp_trialdata.drug]' + 1;
                
                groups = {condition_text(factor.condition)', drug_text(factor.drug)'};
        end
        
        % ANOVA
        [P(iunit,:),statstable{iunit},stats{iunit},terms{iunit}]  = anovan(tmp_rate,groups,'full',3,group_names,'off');
        
    end
    p_anova = array2table(P, 'VariableNames', statstable{iunit}(2:(end-2),1));
        
    savefilename = fullfile(path_target, sprintf('rate_ANOVA_%s.mat', timewin_fields{itw}));
    save(savefilename, 'p_anova', 'statstable', 'stats', 'terms', 'time_windows');
    
end