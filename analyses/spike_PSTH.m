function [ PSTH ] = spike_PSTH( mu, timeEpoch, chann, binSize )
% peristimulus time histogram

    
    numTrial = size(mu,2);
    numBin = ceil( (timeEpoch(2)-timeEpoch(1))/binSize );
    
    % spike number in a particular jitter bin on a particular trial
    PSTH = zeros(numBin, 1);
    
    % only include specified trials
    for iTrial = 1:numTrial
        
        spikes = mu{chann,iTrial, : };
        spikes = spikes( spikes>timeEpoch(1) & spikes<timeEpoch(2));

        bin = ceil( (spikes - timeEpoch(1))/binSize );
  
        PSTH_trial = zeros(numBin,1);
        % estimation based on arbitrary number of spikes per bin
        binT = tabulate(bin);
       
        if (~isempty(binT))
            PSTH_trial(binT(:,1)) = binT(:,2);
        end
        PSTH = PSTH + PSTH_trial;
    end

end

