function RT_drug_modulation(recinfo, trialdata, path_target)
% RT_drug_modulation(recinfo, trialdata, path_target)
%
% Computes average RT across conditions
%
% Parameters
% ----------
% recinfo : table
%     table with single row (one recording)
% trialdata : struct
%     struct with metadata about the trials
% path_target : string
%     path where to store results
% 

% check input 
assert(height(recinfo)==1, 'height recInfo ~= 1')

% remove excluded trials
trialdata([trialdata.exclude]) = [];

% get metadata
num_trial = length(trialdata);
conditions = unique([trialdata.cond_num]);
idx_cond = [trialdata.cond_num]';
idx_drug = [trialdata.drug]' + 1;

% define path to store result
path_target = fullfile(path_target, 'RT_drug_modulation', recinfo.Subject, recinfo.Date);
if ~isfolder(path_target)
    mkdir(path_target)
end

% find normalisation factor
RT_list = [trialdata.RT_EPP];

% store average for each condition
RT = NaN(1,2,length(conditions));
for idrug = 1:2
    for icond = 1:length(conditions)
        
        % get trial indices
        trial_index = (...
            idx_drug==idrug) ...
            & idx_cond==icond;
        
        % pupil
        [RT(1,idrug,icond)] = mean(RT_list(trial_index));
    end
end

savefilename = fullfile(path_target, 'RT.mat');
save(savefilename, 'RT');




