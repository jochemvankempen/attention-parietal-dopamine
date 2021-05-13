function [gain,nlls,nlls2,ps]=data_to_Goris_model(x)
[ps, llike, pvari]= fit_Goris_model(x);

Nc = size(x,1);
mus = zeros(Nc,1);
for j = 1:Nc
    y = x(j,x(j,:)>-1);
    mus(j)= mean(y);
end
% [nlls] = f_loglike_negbin_Goris(x,1/ps,mus);
% [nlls2] = f_loglike_poiss_Goris(x,mus);
[nlls] = f_neg_loglike_negbin_Goris(x,1/ps,mus); %%calculate the negative
%binomial log-likelihood given the estimated gains
[nlls2] = f_neg_loglike_poiss_Goris(x,mus); %%calculate the Posson log-likelihood

gain=ps;
%  if ~isempty(find(nlls<0)) || ~isempty(find(nlls2<0))
%      
%      loglog(nlls2,nlls,'*')
%      hold on
%      loglog([min(nlls2) max(nlls2)],[min(nlls2) max(nlls2)],'k')
% end
