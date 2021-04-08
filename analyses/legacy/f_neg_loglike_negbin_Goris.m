function [nlls] = f_neg_loglike_negbin_Goris(x,r,mus)

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
    nlls(j) = -sum(gammaln(r+y)) + n*gammaln(r) - n*r*log(r/(xbar+r)) - sumx*log(xbar/(xbar+r))+ sum(log(gamma(y+1)));
end

