function [PSTH, time] = spike_convolute_gaussian(mu, time_epoch, chann, sigma, mfactor)
% [PSTH, time] = spike_convolute_gaussian(mu, time_epoch, chann, sigma)
% Convolve spike times with Gaussian, with smooting sigma
%
% Parameters
% ----------
% mu : cell array
%     cell array of size (num_channel, num_trial) that contains spike times
% time_epoch : array
%     time epoch that should be used for the analysis [tStart, tEnd]
% chann : double
%     double indicating channel index
% sigma : float
%     smoothing factor
% mfactor : float
%     multiplication factor. Spike times need to be converted to
%     milliseconds. So for spike times in seconds, mfactor is 1e3. Result
%     is in original units
%
% Returns
% -------
% PSTH : array
%     array of size (num_trial, time) with convolved spikes
% time : array
%     vector of times
%
%


if nargin < 5
    mfactor = 1;
end

num_trial = size(mu,2);

% get time vector
time_epoch = time_epoch*mfactor;
time = time_epoch(1):1:time_epoch(2);

% spike number in a particular bin on a particular trial
PSTH = zeros(num_trial, length(time));

% only include specified trials
for itrial = 1:num_trial

    spikes = mu{chann,itrial};

    PSTH(itrial,:) = convolute_gaussian(spikes*mfactor, time*mfactor, sigma*mfactor);

end



function trialdata = convolute_gaussian(spike_times, time, sigma)
% [trialdata] = spike_convoluteGaussian(spike_times, time, sigma)
% Convolve spike times with Gaussian, with smooting sigma
%
% Example
% -------
% ::
% 
%     [trialdata] = spike_convoluteGaussian(spike_times, 0:0.001:1, 0.01)
%     plot(time, trialdata), plots the histogram for the current trial
%
% Parameters
% ----------
% spike_times : array
%     vector with times at which spike occured
% time : array
%     vector of timepoints
% sigma : float
%     smoothing factor
%
% Returns
% -------
% trialdata : array
%     convolved spikes
%

% convert to row vector
spike_times = spike_times(:)';

trialdata=zeros(1,length(time));
y=ones(1,length(time));
for ispike = 1:length(spike_times)
    spike=(y*(1/(sigma*sqrt(2*pi))).*exp(-(((time-spike_times(ispike)).^2)/(2*sigma^2))));
    trialdata(1,:)=trialdata(1,:)+spike;
end