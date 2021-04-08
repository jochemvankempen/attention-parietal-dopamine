function [parmhat, llike] = nbinfit2(ys,mus,s2s,nts)
%%%%ys: vector with all trials from all conditions
%%%mus: mean for each condition
%%%%s2s: variances for each condition
%%%nts: cumulative number of trials per condition
%%%%parmhat: estimated gain parameter
%%%%llike: negative log-likelihood for the estimated parameter
Nc = length(mus);
options = statset('nbinfit');
if ~isfloat(ys)
    ys = double(ys);
end

% % Ensure that a NB fit is appropriate.
% for i = 1:Nc
%     xbar = mus(i);
%     s2 = s2s(i);
%     if s2 <= xbar
%         parmhat = cast([Inf 1.0],class(ys));
%         parmci = cast([Inf 1; Inf 1],class(ys));
%         warning('stats:nbinfit:MeanExceedsVariance',...
%             'The sample mean exceeds the sample variance -- use POISSFIT instead.');
%         return
%     end
% end


% Use Method of Moments estimates as starting point for MLEs.
rhats = zeros(Nc,1);
for i = 1:Nc
    xbar = mus(i);
    s2 = s2s(i);
    rhats(i) = (xbar.*xbar) ./ (s2-xbar);
end
rhat = nanmean(rhats);

if rhat < 0
    parmhat = cast([NaN],class(ys));
    llike = cast([NaN],class(ys));
    warning('stats:nbinfit:MeanExceedsVariance',...
        'The sample mean exceeds the sample variance -- use POISSFIT instead.');
    return
end

% Parameterizing with mu=r(1-p)/p makes this a 1-D search for rhat.
[rhat,nll,err,output] = ...
    fminsearch(@negloglike, rhat, options, nts, ys, mus, options.TolBnd);
if (err == 0)
    % fminsearch may print its own output text; in any case give something
    % more statistical here, controllable via warning IDs.
    if output.funcCount >= options.MaxFunEvals
        wmsg = 'Maximum likelihood estimation did not converge.  Function evaluation limit exceeded.';
    else
        wmsg = 'Maximum likelihood estimation did not converge.  Iteration limit exceeded.';
    end
    if rhat > 100 % shape became very large
       wmsg = sprintf('%s\n%s', wmsg, ...
                      'The Poisson distribution might provide a better fit.');
    end
    warning('stats:nbinfit:IterOrEvalLimit',wmsg);
elseif (err < 0)
    error('stats:nbinfit:NoSolution', ...
          'Unable to reach a maximum likelihood solution.');
end
parmhat = [1/rhat];
llike = nll;


%-------------------------------------------------------------------------

function nll = negloglike(r, nts, y, mus, tolBnd)
% Objective function for fminsearch().  Returns the negative of the
% (profile) log-likelihood for the negative binomial, evaluated at
% r.  From the likelihood equations, phat = rhat/(xbar+rhat), and so the
% 2-D search for [rhat phat] reduces to a 1-D search for rhat -- also
% equivalent to reparametrizing in terms of mu=r(1-p)/p, where muhat=xbar
% can be found explicitly.

% Restrict r to the open interval (0, Inf).
Nc = length(mus);
nlls = zeros(Nc,1);
if r < tolBnd
    nll = Inf;
else
    for i=1:Nc
      x = y(1+nts(i):nts(i+1));  
      sumx = sum(x);
      n = length(x);
      xbar = mean(x);
      nlls(i) = -sum(gammaln(r+x)) + n*gammaln(r) - n*r*log(r/(xbar+r)) - sumx*log(xbar/(xbar+r));
    end
    nll = sum(nlls);
end







