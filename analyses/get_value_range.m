function x_pos = get_value_range(x, proportion)
% x_ratio = get_value_range(x, proportion)
%
% Get the value at position proportion in the range spanned by x. 
%
% Parameters
% ----------
% x : array
%     array of numeric values
% proportion : float
%     float indicating the position for which the value on range(x) should be
%     returned. This is value between zero and 1.
% 
% Returns
% -------
% x_ratio : float
%     float indicating the value of the position for which the value on range(x) should be
%     returned
%
% Example
% -------
% ::
% 
%     x_pos = get_value_range([-3 10], 0.6)
%     x_pos =
%         4.8000
% 

assert((0<=proportion) & (proportion<=1), 'proportion should be between zero and 1')

x_pos = (max(x)-min(x))*proportion+min(x);
