function [ps, llike, pvari] = fit_Goris_model(x)

%%%%x: input matrix with conditions(rows) and trials(columns); negative values
%%%%  are used to indicate different number of trials per conditions
%%%%ps: parameter of gain of the gamma distribution. A single parameter is 
%%%%inferred and assumed constant for all conditions
%%%%llike: total negative log-likelihood for the selected parameter value
%%%%pvari: percentage of variance (normalized to one) explained by the gain fluctuations for
%%%%each condition. There are two columns, the first follows eq. 3 of Goris
%%%%et al 2014. The second normalizes by the real variance.

Nc = size(x,1);
mus = zeros(Nc,1);%%means for each condition
s2s = zeros(Nc,1);%%%variance for each condition
nts = zeros(Nc,1);%%%number of trials for each condition
ys = [];
for j = 1:Nc
  y = x(j,x(j,:)>-1);  
  mus(j) = mean(y);
  s2s(j) = var(y);
  ys = [ys y];
  nts(j)=length(y);
end
nts2 = [0 cumsum(nts)'];
% figure
% plot(mus,s2s,'*')
% hold on
% plot([0:1:max(mus)], [0:1:max(mus)],'k')
% xlabel('mean')
% ylabel('variance')
pvari = zeros(length(mus),2);
[ps, llike] = nbinfit2(ys,mus,s2s,nts2); %%%fitting of the negative binomial.
pvari(:,1) = (ps*mus.^2)./(mus+ps*mus.^2);
pvari(:,2) = (ps*mus.^2)./(s2s);


