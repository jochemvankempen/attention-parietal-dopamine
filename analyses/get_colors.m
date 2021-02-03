function colors = get_colors(type)
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

        colors(1,1,:) = [0 0 102]/255 ;
        colors(1,2,:) = [0 128 255]/255 ;
        colors(2,1,:) = [153 0 0]/255 ;
        colors(2,2,:) = [255 128 0]/255 ;
        
    case 'drug'

        colors = [...
            0.2 0.2 0.2;
            0.5 0.5 0.5
            ];
end



