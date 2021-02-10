function [xspikes, yspikes] = spike_raster( mu, timeEpoch, chann, yoffset )
%
% create matrix with raster
%
% Example
% -------
% ::
% 
%     [xspikes, yspikes] = spike_compute_raster(sps, [-0.1 0.3], 2, 10)
%     plot(xpsikes, yspikes)
% 
% Parameters
% ----------
% sps : cell
%     cell array of size (numChannel x numTrial) with spike times relative
%     to event
% timeEpoch : array
%     array of size (1x2) with start and endpoint of time epoch.
% chann : double
%     double indicating channel number
% yoffset : double
%     scalar that indicates the offset of the raster (y-axis), default is 0
%
% Returns
% -------
% xspikes : array 
%     array of size (numTrial
% yspikes : array
% 
% 
% 

if nargin<4
    yoffset=0;
end

[numChann, numTrial] = size(mu);


xspikes = [];
yspikes = [];

for itrial = 1:numTrial
    sps = mu{chann, itrial};
    
    if ~isempty(timeEpoch)
        sps = sps( sps>timeEpoch(1) & sps<timeEpoch(2) );
    end
    
    [xtmp, ytmp] = compute_raster(sps, (itrial + yoffset) );
    
    xspikes = [xspikes xtmp];
    yspikes = [yspikes ytmp];
end

function [xspikes, yspikes] = compute_raster(spike_times, y, y_length)
% computes a raster from spike times
% 
% INPUT:
% spike_times: vector with times at which spike occured
% y: y-position, where to plot this trials' spikes (e.g. trial/channel index)
% y_length: The length of the spike to be drawn
% 
% OUTPUT:
% xspikes, yspikes: x and y values of the spikes
%
% plot(xspikes, yspikes) plots the raster for the current trial
%
% Jochem van Kempen, 2017
 
if nargin<3
    y_length = 1;
end

% convert to row vector
spike_times = spike_times(:)';

xspikes = repmat(spike_times,3,1);
yspikes = nan(size(xspikes));
if ~isempty(yspikes)
    yspikes(1,:) = y-y_length;
    yspikes(2,:) = y;
end