function [stats2report, predict_data, STATS, model] = fitlme_singleVar_sequential(data_table, plot_fit, onlyLinearFit)
% sequentially fit LME to data and compare whether linear or quadratic fit is better
% than constant or linear fit, respectively.
% 
% Parameters
% ----------
% data_table : table 
%     table with variables 'y' and 'x' as dependent and independent
%     variable
% plot_fit : boolean
%     plot the fit of the model, default is 0
% onlyLinearFit : boolean
%     only apply linear fit, default is 0
%
% Returns
% -------
% stats2report : table
%     table with (LRStat) loglikelihood of model comparison, (delta DF) the
%     difference in degrees of freedom in the models that were compared and
%     (pValue) the p value of model comparison 
% predict_data : struct
%     struct with the mean fit of the model, CIFit: the confidence intervals of the model fit
% CIFit: the standard error intervals of the model fit
% STATS: stats on the model, the raw beta values etc. 
%
% Example:
% datamat = table(RT, Subject, Bin, 'VariableNames', {'Y','Subject','Bin'});
% [stats2report, meanFit, CIFit, SEFit, STATS, model] = fitlme_singleVar_sequential(datamat, 0);
%
% Jochem van Kempen, 26/06/2018

if nargin < 2
    plot_fit = 0;
end

if ~exist('onlyLinearFit','var')
    onlyLinearFit = 0;
end

lme = fitlme(data_table,['y ~ 1 '], 'FitMethod', 'ML');% fit constant
lme2 = fitlme(data_table,['y ~ x'], 'FitMethod', 'ML'); % fit linear model
if ~onlyLinearFit
    lme3 = fitlme(data_table,['y ~ x^2'], 'FitMethod', 'ML'); % fit quadratic model (this includes the linear model), note this doesn't have orthogonal contrasts!!
end
% compare the model fits
comp1 = compare(lme, lme2, 'CheckNesting',true); % compare linear to constant
if ~onlyLinearFit
    comp2 = compare(lme2, lme3, 'CheckNesting',true); % compare quadratic to linear
else
    comp2.pValue = 1;
end

x_range = [min(data_table.x), max(data_table.x)];
x_range = x_range + round([-mean(x_range)/2 mean(x_range)/2]);
predict_data = table((x_range(1):x_range(2))', 'VariableNames', {'x'});

if comp2.pValue < 0.05
    [~,~,STATS] = fixedEffects(lme3);    
    [fit, fit_CI] = predict(lme3, predict_data); % get the fitted values for the new x values
    stats2report = array2table([comp2.LRStat(2) comp2.deltaDF(2) comp2.pValue(2)], 'VariableNames', {'LRStat','deltaDF','pValue'});
    STATS = STATS(3,:);
    model = 'quadratic';
        
elseif comp1.pValue < 0.05
    [fit, fit_CI] = predict(lme2, predict_data); % get the fitted values for each subject
    [~,~,STATS] = fixedEffects(lme2);
    stats2report = array2table([comp1.LRStat(2) comp1.deltaDF(2) comp1.pValue(2)], 'VariableNames', {'LRStat','deltaDF','pValue'});
    STATS = STATS(2,:);
    model = 'linear';
    
else
    stats2report = array2table([comp1.LRStat(2) comp1.deltaDF(2) comp1.pValue(2)], 'VariableNames', {'LRStat','deltaDF','pValue'});
    fit_mean = [];
    fit_CI = []; fit_SE = []; STATS = [];
    model = 'none';
    return
    
end

% fit_Y = fit;%reshape(fit, [nSubject,nBin]);% reshape back into subject x Bin
% fit_Y_CI1 = fitCI(:,1);%reshape(fitCI(:,1), [nSubject,nBin]);% reshape back into subject x Bin
% fit_Y_CI2 = fitCI(:,2);%reshape(fitCI(:,2), [nSubject,nBin]);% reshape back into subject x Bin

% fit_mean = fit;
% fit_CI = [1:length(fit) length(fit):-1:1 ; fit_CI(:,1)' flip(fit_CI(:,2))'];
% fit_SE = [];
% SEFit = [1:nBin nBin:-1:1; meanFit - std(fit_Y)/sqrt(nSubject) flip(meanFit + std(fit_Y)/sqrt(nSubject))];

predict_data.fit_mean = fit;
predict_data.fit_CI_low = fit_CI(:,1);
predict_data.fit_CI_high = fit_CI(:,2);


if plot_fit
    var2plot = reshape(eval(['dataMat.' y ]), [nSubject,nBin]);
    figure
    hold on
    h = plot(fit_mean, 'linewidth', 2);
    h = patch(fit_CI(1,:), fit_CI(2,:), h.Color);
    h.FaceAlpha = 0.3;
    h.EdgeAlpha = 0.3;
    h.EdgeColor = [h.FaceColor];
    plot(mean(var2plot),'.','markersize',15)
end






