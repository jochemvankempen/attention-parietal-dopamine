function [colors, labels, labels2] = get_colors(type)
% colors = get_colors(type)
%
% Get the colors used for plotting
% 
% Parameters
% ----------
% type : string
%     string indicating which colors to use
%
% Returns
% -------
% colors : array of floats
%     colors to plot, color is last dimension
%

switch type
    
    case 'att_drug'

%         % colour
        colors(1,1,:) = [0 0 102]/255 ;
        colors(1,2,:) = [0 128 255]/255 ;
        colors(2,1,:) = [153 0 0]/255 ;
        colors(2,2,:) = [255 128 0]/255 ;
        
%         % grayscale
%         colors(1,1,:) = [50 50 50]/255 ;
%         colors(1,2,:) = [50 50 50]/255 ;
%         colors(2,1,:) = [153 153 153]/255 ;
%         colors(2,2,:) = [153 153 153]/255 ;
        
        labels{1,1} = 'Attend RF, Drug off';
        labels{1,2} = 'Attend RF, Drug on';
        labels{2,1} = 'Attend away, Drug off';
        labels{2,2} = 'Attend away, Drug on';
        
        labels2{1,1} = 'RF, no drug';
        labels2{1,2} = 'RF, drug';
        labels2{2,1} = 'Away, no drug';
        labels2{2,2} = 'Away, drug';
        
    case 'drug'

        colors = [...
            0.1 0.1 0.1;
            0.5 0.5 0.5
            ];
        
    case 'spikewidth'
%         colors(1,:) = [0 0 0]/255 ;
%         colors(2,:) = [153 153 153]/255 ;

        colors(1,:) = [0        0.4470  0.7410];
        colors(2,:) = [0.8500   0.3250  0.0980];
end



