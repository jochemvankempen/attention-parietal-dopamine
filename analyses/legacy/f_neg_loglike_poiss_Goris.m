function [nlls] = f_neg_loglike_poiss_Goris(x,mus)

%%%%x input matrix with conditions x trials
%%%%negative values are used to indicate different number of trials per
%%%%conditions
Nc = size(x,1);
nlls = zeros(Nc,1);

for j = 1:Nc
    y = x(j,x(j,:)>-1);
    sumx = sum(y);
    n = length(y);
    xbar = mus(j);
    nlls(j) = + n*xbar - sumx*log(xbar)+ sum(log(gamma(y+1)));
end