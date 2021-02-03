function [ area, hit, fa ] = ROC_area(data, class)
% [ area, hit, fa ] = ROC_area(data, class )
%
% Computes the area under the ROC curves
%
% Parameters
% ----------
% data : array of float
%     vector containing a set of samples
% class : array of boolean
%     vector with a membership of each data point in the class 0 or 1 (binary only)
% 
% Returns
% -------
% area : array of bool
%     booleans indicating trial inclusion
% 

% perform checks
assert(length(data)~=length(class), 'Data and their class labels should have the same length');
assert(all((class==0 | class==1)), 'Class labels should be either 0 or 1');

n = length(data);

[dataSort, indSort] = sort(data);
classSort = class(indSort);

hit = zeros(n,1);
fa = zeros(n,1);

n1 = sum( classSort==1 );
n2 = sum( classSort==0 );
% for each value of the threshold
for iThresh = 1:n
    hit(n-iThresh+1) = sum( classSort(iThresh:n)==1 )/n1;  % number of hits
    fa(n-iThresh+1)  = sum( classSort(iThresh:n)==0 )/n2;   % number of fals alarms
end

deltaFA = fa(2:end)-fa(1:end-1);
area = deltaFA'*hit(2:end);

end

