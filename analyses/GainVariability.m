classdef GainVariability
    % Implementation of gain variability estimation using a negative
    % binomial distribution fit to spike count/rate, according to Goris et
    % al., 2014 Nature Neuroscience doi:﻿10.1038/nn.3711.
    
    properties
              
        num_cond; % number of conditions
        num_trial; % number of trials for each condition
        num_trial_cum; % cumulative number of trials for each condition
        mus; % means for each condition
        s2s; % variance for each condition
        spike_concat; % concatenated spike count/rate across conditions
        
        rhat; % gain parameter of the gamma distribution. This is the inverse of the rate parameter of the negative binomial distribution. 
        negloglike;
        pvariance;
    end
    
    
    methods
        
        function obj = GainVariability(X)
            % obj = GainVariability(X)
            %
            % class constructor
            %
            % Parameters
            % ----------
            % X : array
            %   array of size (conditions x trials) with spike count/rate.
            %   Negative values are used to indicate different number of
            %   trials per condition
            % 
            
            % number of conditions
            [obj.num_cond, tmp_trial] = size(X);
            assert(tmp_trial>obj.num_cond, 'Number of trials is smaller than number of conditions, check input dimensions of X')
            
            % init
            obj.mus = zeros(obj.num_cond,1); % means for each condition
            obj.s2s = zeros(obj.num_cond,1); % variance for each condition
            obj.num_trial = zeros(obj.num_cond,1); % number of trials for each condition
            
            % compute mu, variance and number of trials for each condition
            obj.spike_concat = [];
            for j = 1:obj.num_cond
                spike_cond = X(j,X(j,:)>-1); % get spike count/rate for each condition, omit trials with negative values
                obj.mus(j) = mean(spike_cond);
                obj.s2s(j) = var(spike_cond);
                obj.spike_concat = [obj.spike_concat spike_cond];
                obj.num_trial(j)=length(spike_cond);
            end
            obj.num_trial_cum = [0 cumsum(obj.num_trial)'];
            
            % figure
            % plot(obj.mus,obj.s2s,'*')
            % hold on
            % plot([0:1:max(obj.mus)], [0:1:max(obj.mus)],'k')
            % xlabel('mean')
            % ylabel('variance')
        end
        
        function obj = fit(obj)
            % obj = fit(obj)
            % 
            % Fit the negative binomial distribution to spike count/rate to
            % acquire the gain variability parameter. A single parameter is
            % inferred and assumed constant for all conditions. 
            % 
            % Returns
            % -------
            % paramhat : float
            %   gain parameter of the gamma distribution. This is the
            %   inverse of the rate parameter of the negative binomial
            %   distribution.
            % nll : float
            %   float indicating the total negative log-likelihood for the
            %   selected parameter value 
            % pvariance : array of floats
            %   percentage of variance (normalized to one) explained by the
            %   gain fluctuations for each condition. There are two
            %   columns, the first follows eq. 3 of Goris et al 2014. The
            %   second normalizes by the real variance. 
            %
            
            % fitting of the negative binomial.
            obj = nbinfit(obj); 
            
            % computing the percentage variance explained
            obj.pvariance = zeros(length(obj.mus),2);
            obj.pvariance(:,1) = (obj.rhat*obj.mus.^2)./(obj.mus+obj.rhat*obj.mus.^2); % percentage of variance according to eq3 of Goris et al.
            obj.pvariance(:,2) = (obj.rhat*obj.mus.^2)./(obj.s2s); % percentage of variance normalised by the real variance.
            
        end
        
        
        function obj = nbinfit(obj)
            
            %%%%ys: vector with all trials from all conditions
            %%%mus: mean for each condition
            %%%%s2s: variances for each condition
            %%%nts: cumulative number of trials per condition
            %%%%parmhat: estimated gain parameter
            %%%%llike: negative log-likelihood for the estimated parameter

            options = statset('nbinfit');
            if ~isfloat(obj.spike_concat)
                obj.spike_concat = double(obj.spike_concat);
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
            paramhats = zeros(obj.num_cond,1);
            for icond = 1:obj.num_cond
                xbar = obj.mus(icond);
                s2 = obj.s2s(icond);
                paramhats(icond) = (xbar.*xbar) ./ (s2-xbar);
            end
            paramhat = nanmean(paramhats);
            
            if paramhat < 0
                obj.rhat = cast([NaN],class(obj.spike_concat));
                obj.negloglike = cast([NaN],class(obj.spike_concat));
                warning('stats:nbinfit:MeanExceedsVariance',...
                    'The sample mean exceeds the sample variance -- use POISSFIT instead.');
                return
            end
            
            % Parameterizing with mu=r(1-p)/p makes this a 1-D search for rhat.
            [paramhat,nll,err,output] = ...
                fminsearch(@negloglike_nbin, paramhat, options, obj, options.TolBnd);
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
            
            % store
            obj.rhat = 1/paramhat;
            obj.negloglike = nll;            
        
        end
        
        
        function [nll,nlls] = negloglike_nbin(r, obj, tolBnd)
            % Objective function for fminsearch().  Returns the negative of the
            % (profile) log-likelihood for the negative binomial, evaluated at
            % r.  From the likelihood equations, phat = rhat/(xbar+rhat), and so the
            % 2-D search for [rhat phat] reduces to a 1-D search for rhat -- also
            % equivalent to reparametrizing in terms of mu=r(1-p)/p, where muhat=xbar
            % can be found explicitly.
            
            % init
            nlls = zeros(obj.num_cond,1);
            
            % Restrict r to the open interval (0, Inf).
            if r < tolBnd
                nll = Inf;
            else
                for icond=1:obj.num_cond
                    
                    % extract trials of this condition
                    x = obj.spike_concat(1+obj.num_trial_cum(icond):obj.num_trial_cum(icond+1));
                    
                    % get sum, length and mean of x
                    sumx = sum(x);
                    n = length(x);
                    xbar = mean(x);
                    
                    % negative log-likelihood per condition
                    nlls(icond) = -sum(gammaln(r+x)) + n*gammaln(r) - n*r*log(r/(xbar+r)) - sumx*log(xbar/(xbar+r));
                end
                
                % negative log-likelihood across all conditions
                nll = sum(nlls);
                
            end
        end
        
        function [nll,nlls] = negloglike_poiss(X,mus)
            % [nll,nlls] = negloglike_poiss(x,mus)
            % 
            % compute negative log-likelihood for each condition (row) of
            % data given Poisson distribution with param: mus (lambda)
            % 
            % Parameters
            % ----------
            % X : array
            %   array of size (conditions x trials) with spike count/rate.
            %   Negative values are used to indicate different number of
            %   trials per condition
            % mus : array
            %   array of length (conditions) with mean spike count/rate for
            %   each condition
            %
            % Returns
            % -------
            % nll : float
            %   float indicating the negative log-likelihood across
            %   conditions
            % nlls : array of floats
            %   floats indicating the negative log-likelihood for each
            %   condition
            %

            % number of conditions
            num_cond = size(X,1);
            assert(num_cond==length(mus), 'inconsistent condition numbers')
            
            % init
            nlls = zeros(num_cond,1);
            
            % loop over conditions and compute negative log-likelihood
            for icond = 1:num_cond
                
                % extract trials of this condition, omit trials with
                % negative spike rate
                x = X(icond,X(icond,:)>-1);
                
                % get sum, length and mean of x
                sumx = sum(x);
                n = length(x);
                xbar = mus(icond);
                
                % negative log-likelihood per condition
                nlls(icond) = n*xbar - sumx*log(xbar)+ sum(log(gamma(x+1)));
            end
            
            % negative log-likelihood across all conditions
            nll = sum(nlls);

        end
        
    
    end
end

