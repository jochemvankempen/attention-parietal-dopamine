function [P, STATS, textOut, fullTextOut] = compare_means(Mat1, Mat2, paired, forceMethod)
% [P, STATS, textOut, fullTextOut] = compare_means(Mat1, Mat2, paired, forceMethod)
% compute stats to compare mean values to each other
%
% Parameters
% ----------
% Mat1 : array
%     array of size (n, 1) with values to compare to Mat2
% Mat2 : array
%     array of size Mat1 or single value to compare Mat1 to (default, 0)
% paired : boolean 
%     boolean indicating whether to perform a paired (1) or unpaired (0,
%     default) test. 
% forceMethod : string (optional)
%     'ttest' or 'rank', to force usage of certain statistical test
%
%
if nargin<4 || isempty(forceMethod)
    forceMethod = 0;
else
    switch forceMethod
        case 'ttest'
            forceMethod = 1;
        case 'rank'
            forceMethod = 2;
    end
end

if nargin < 3 || isempty(paired)
    paired = 0;
end

if nargin < 2 || isempty(Mat2)
    Mat2    = 0;
end

if length(Mat2) == 1
    paired  = 1;
end

% test for normality
Hnorm = kstest(Mat1);

Hnorm2 = 0;
if length(Mat2) > 1
    Hnorm2 = kstest(Mat2);
end

if (forceMethod==2) || (Hnorm || Hnorm2)
    % if non-normal, 
    if paired
        [P, H, STATS] = signrank(Mat1, Mat2, 'method', 'exact');
    else
        [P, H, STATS] = ranksum(Mat1, Mat2, 'method', 'exact'); 
    end
    fullTextOut = sprintf('ws(%1.1f), \n  p = %1.1e', STATS.signedrank, P);
else

    if paired
        [H, P, CI, STATS] = ttest(Mat1, Mat2);
    else
        [H, P, CI, STATS] = ttest2(Mat1, Mat2);
    end
    fullTextOut = sprintf('t(%d) = %1.1f, \n  p = %1.1e', STATS.df, STATS.tstat, P);
    
end

if P<0.001
    %textOut = sprintf('p = %1.1e', P);
    b = floor(log10(abs(P))) + 1;
    textOut = sprintf('p < 10^{%d}', b);
else
    textOut = sprintf('p = %1.3f', P);
end









