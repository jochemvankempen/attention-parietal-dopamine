function spike_rate_summary(recinfo, trialdata, unit, time_windows, PSTH_windows, path_target)
% spike_rate_summary(recinfo, spike_rate, path_target)
%
% Computes average rate and fano factor for each condition over each field in time_windows
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
% PSTH_windows : struct
%     struct with time windows for which to extract PSTH
% path_target : string
%     path where to store results
% 

% check input 
assert(height(recinfo)==1, 'height recInfo ~= 1')

% get dimensions
[num_unit, num_trials] = size(unit.StimAlign);
conditions = unique([trialdata.cond_num]);
num_group = length(unique([trialdata.group]));

% compute firing rate
rate = spike_rate(unit, time_windows, 'rate');

% get alignments
timewin_fields = fields(rate);
PSTH_fields = fields(PSTH_windows);

% define path to store result
path_target = fullfile(path_target, 'spike_rate_summary', recinfo.Subject, recinfo.Date);
if ~isfolder(path_target)
    mkdir(path_target)
end

% compute histogram
binsize = 30;
maxhist = zeros(num_unit,1); % normalisation factor
for itw = 1:length(PSTH_fields)
        
    PSTH_event = PSTH_windows.(PSTH_fields{itw}){1};
    time_window = PSTH_windows.(PSTH_fields{itw}){2};
    
    for iunit = 1:num_unit
        
        % select only trial window for which this unit has spikes
        [tmp_trialdata, tmp_unit] = remove_excluded_trials(trialdata, unit, iunit);

        % compute histogram
        [hist_spike, hist_time] = spike_convolute_gaussian(tmp_unit.([PSTH_event 'Align']), time_window, iunit, binsize);
        
        % init store
        if iunit==1
            PSTH.(PSTH_fields{itw}).samples = NaN(num_unit, 2, length(conditions), length(hist_time)); % unit, drug, conditions, time
            PSTH.(PSTH_fields{itw}).time = hist_time;
        end
        
        % rate across condition
        idx_cond = [tmp_trialdata.cond_num]';
        idx_drug = [tmp_trialdata.drug]' + 1;
        
        for idrug = 1:2
            for icond = 1:length(conditions)
                
                % get trial indices
                idx_trial = (...
                    idx_drug==idrug) ...
                    & idx_cond==icond;
                
                % store
                PSTH.(PSTH_fields{itw}).samples(iunit,idrug,icond,:) = mean(hist_spike(idx_trial,:),1);
                
                if maxhist(iunit)<max(PSTH.(PSTH_fields{itw}).samples(iunit,idrug,icond,:))
                    maxhist(iunit) = max(PSTH.(PSTH_fields{itw}).samples(iunit,idrug,icond,:));
                end
            end
        end
    end
end

% normalisation
for itw = 1:length(PSTH_fields)
    for iunit = 1:num_unit
        PSTH.(PSTH_fields{itw}).samples(iunit,:,:,:) = PSTH.(PSTH_fields{itw}).samples(iunit,:,:,:) / maxhist(iunit) * 100; 
    end
end
savefilename = fullfile(path_target, sprintf('rate_PSTH.mat'));
save(savefilename, 'PSTH', 'PSTH_windows', 'maxhist');


% loop over time windows
for itw = 1:length(timewin_fields)   
    
    % init
    [rate_cond, FF_cond] = deal(NaN(num_unit, 2, length(conditions), num_group));% unit, drug, cond, group
    [lm.coef] = deal(NaN(num_unit, 2, 2)); % unit, drug, coefficients
    [lm.tStat, lm.pValue] = deal(NaN(num_unit, 2)); % unit, drug
    
    for iunit = 1:num_unit
                
        % select only trial window for which this unit has spikes
        [tmp_trialdata, tmp_unit, idx_exclude_trials] = remove_excluded_trials(trialdata, unit, iunit);
        tmp_rate = rate.(timewin_fields{itw})(iunit,~idx_exclude_trials)';
        
        % get trial indices
        idx_cond = [tmp_trialdata.cond_num]';
        idx_drug = [tmp_trialdata.drug]' + 1;
        idx_group = [tmp_trialdata.group]';
        
        block_change = [tmp_trialdata.block_change]';

        for idrug = 1:2
            
            % rate across condition
            for icond = 1:length(conditions)
                for igroup = 1:num_group
                    
                    % get trial indices
                    idx_trial = ...
                        idx_drug==idrug ...
                        & idx_cond==icond ...
                        & idx_group==igroup;
                    
                    % rate
                    [rate_cond(iunit,idrug,icond,igroup)] = mean(tmp_rate(idx_trial));
                    
                    % FF
                    [FF_cond(iunit,idrug,icond,igroup)] = var(tmp_rate(idx_trial)) / mean(tmp_rate(idx_trial));
                    
                end
            end
            
        end
        
%         rate_normalisation = norm(tmp_rate);
        rate_normalisation = 1;
        for idrug = 1:2
            % rate over time/blocks
            
            % exclude first block, as no drug has been applied yet
            idx_block = true(length(tmp_trialdata),1);
            idx_block_exclude = find(diff([tmp_trialdata.block_num]),1,'first');
            idx_block(1:idx_block_exclude) = false;
            
            % trial indices
            idx_trial = ...
                idx_drug==idrug ...
                & idx_block;

            % fit lm
            rate_normalised = tmp_rate(idx_trial)/rate_normalisation;
            t = table(block_change(idx_trial), rate_normalised, 'VariableNames', {'x','y'});
            try
                mfit = fitlm(t, 'linear');
                % store results
                lm.coef(iunit,idrug,:) = mfit.Coefficients.Estimate;
                lm.tStat(iunit,idrug) = mfit.Coefficients.tStat(2);
                lm.pValue(iunit,idrug) = mfit.Coefficients.pValue(2);
            catch
                warning('Could not fit lm to block_change x rate')
            end
            
%             % plot
%             print_text = {
%                 sprintf('Drug: %d', idrug-1)
%                 sprintf('Coef: %1.2e %s %1.2e', mfit.Coefficients.Estimate(2), '\pm', mfit.Coefficients.SE(2)),
%                 sprintf('p = %1.2e', mfit.Coefficients.pValue(2))
%                 };
%             scatter(t.x, t.y)
%             lsline
%             text(25, max(t.y)*0.8, print_text, 'Interpreter', 'tex', 'FontSize', 14)
%             pause
        end
        
        clear tmp*
    end
        
    savefilename = fullfile(path_target, sprintf('rate_summary_%s_%s.mat', timewin_fields{itw}, trialdata(1).group_label));
    save(savefilename, 'rate_*', 'FF_*', 'time_windows', 'lm');
    
end