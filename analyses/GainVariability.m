classdef GainVariability
    % Implementation of gain variability estimation using a negative
    % binomial distribution fit to spike count/rate, according to Goris et
    % al., 2014 Nature Neuroscience doi:﻿10.1038/nn.3711.
    
    properties
        
        num_cond; % number of conditions
        num_trial; % number of trials for each condition
        num_trial_cum; % cumulative number of trials for each condition
        mus; % mean for each condition
        s2s; % variance for each condition
        spike_concat; % concatenated spike count/rate across conditions
        
        rhat; % gain parameter of the gamma distribution. This is the inverse of the rate parameter (dispersion) of the negative binomial distribution.
        negloglike; % negative log-likelihood of negative binomial fit.
        negloglike_cond; % negative log-likelihood of negative binomial fit for each condition.
        negloglike_poiss; % negative log-likelihood of Poisson fit.
        negloglike_poiss_cond; % negative log-likelihood of Poisson fit for each condition.
        pvariance; % proportion variance explained
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
            
            % allow return to initialise empty object
            if nargin==0
                return
            end
            
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
        end
                
        function [m2v_table,fighandle] = eval_mean2variance(obj,xval,showplot)
            % [m2v_table,fighandle] = eval_mean2variance(obj,xval,showplot)
            % 
            % Compute mean to variance relationship with negative binomial
            % fit.
            %
            % Parameters
            % ----------
            % obj : object
            %   GainVariability class object
            % xval : array of int (optional)
            %   array of integers indicating at which x values (mean firing
            %   rates) the variance should be evaluated. Default is
            %   xval=1:100; 
            % showplot : bool (optional)
            %   boolean indicating whether to plot the pmf, default is
            %   false.
            %
            % Returns
            % -------
            % m2v_table : table
            %   table with y (estimated variance of firing rate given
            %   nbinom fit) evaluated at x (mean rate) values. 
            % fighandle : figure handle
            %
            
            % check input.
            if nargin<2 || isempty(xval)
                xval = 1:100;
            end
            if nargin<3
                showplot=false;
            end
            if isempty(obj.rhat)
                obj = obj.fit; % fit negative binomial
            end
            
            % variance of negative binomial distribution
            yval = xval + obj.rhat*xval.^2;
            
            % store in table
            m2v_table = table(xval(:),yval(:),'VariableNames',{'x','y'});
            
            if showplot
                % mean to variance relationship
                figure(1),clf
                hold on
                plot(obj.mus,obj.s2s,'o','MarkerFaceColor','auto')
                h = plot(xval,yval);
                xlim([0.8 max(yval)*1.1])
                ylim([0.8 max(yval)*1.1])
                plot(xlim,ylim,'k','linew',1)
                set(gca,'xscale','log','yscale','log')
                fighandle = gcf;
            else
                fighandle = [];
            end
                        
        end
                
        function obj = fit(obj)
            % obj = fit(obj)
            %
            % Fit the negative binomial distribution to spike count/rate to
            % acquire the gain variability parameter. A single parameter is
            % inferred and assumed constant for all conditions.
            %
            % Parameters
            % ----------
            % obj : object
            %   GainVariability class object
            %
            % Returns
            % -------
            % obj.rhat : float
            %   gain parameter of the gamma distribution. This is the
            %   inverse of the rate parameter of the negative binomial
            %   distribution.
            % obj.negloglike : float
            %   float indicating the total negative log-likelihood for the
            %   selected parameter value
            % obj.pvariance : array of floats
            %   percentage of variance (normalized to one) explained by the
            %   gain fluctuations for each condition. There are two
            %   columns, the first follows eq. 3 of Goris et al 2014. The
            %   second normalizes by the real variance.
            %
            
            % fitting of the negative binomial.
            obj = nbinfit(obj);
            
            % get log likelihood
            options = statset('nbinfit');
            [obj.negloglike,obj.negloglike_cond] = negloglike_nbin(1/obj.rhat,obj,options.TolBnd);
            [obj.negloglike_poiss,obj.negloglike_poiss_cond] = negloglike_poisson(obj);
            
            % computing the percentage variance explained
            obj.pvariance = zeros(length(obj.mus),2);
            obj.pvariance(:,1) = (obj.rhat*obj.mus.^2)./(obj.mus+obj.rhat*obj.mus.^2); % percentage of variance according to eq3 of Goris et al.
            obj.pvariance(:,2) = (obj.rhat*obj.mus.^2)./(obj.s2s); % percentage of variance normalised by the real variance.
            
        end
        
        
        function obj = nbinfit(obj)
            % obj = nbinfit(obj)
            %
            % Called by GainVariability.fit
            % Fit the negative binomial distribution to spike count/rate to
            % acquire the gain variability parameter. A single parameter is
            % inferred and assumed constant for all conditions.
            %
            % Parameters
            % ----------
            % obj : object
            %   GainVariability class object
            %
            % Returns
            % -------
            % obj.rhat : float
            %   gain parameter of the gamma distribution. This is the
            %   inverse of the rate parameter of the negative binomial
            %   distribution.
            % obj.negloglike : float
            %   float indicating the total negative log-likelihood for the
            %   selected parameter value
            %
            % based on built-in MATLAB function nbinfit
            %
            
            options = statset('nbinfit');
            if ~isfloat(obj.spike_concat)
                obj.spike_concat = double(obj.spike_concat);
            end
            
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
                if paramhat > 100 % shape became very large
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
            % [nll,nlls] = negloglike_nbin(r, obj, tolBnd)
            %
            % Objective function for fminsearch().  Returns the negative of the
            % (profile) log-likelihood for the negative binomial, evaluated at
            % r.  From the likelihood equations, phat = rhat/(xbar+rhat), and so the
            % 2-D search for [rhat phat] reduces to a 1-D search for rhat -- also
            % equivalent to reparametrizing in terms of mu=r(1-p)/p, where muhat=xbar
            % can be found explicitly.
            %
            % based on built-in MATLAB function nbinfit/negloglike 
            %
            
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
                    nlls(icond) = -sum(gammaln(r+x)) + n*gammaln(r) - n*r*log(r/(xbar+r)) - sumx*log(xbar/(xbar+r)) + sum(log(gamma(x+1)));
                    %                     nlls(icond) = -sum(gammaln(r+x)) + n*gammaln(r) - n*r*log(r/(xbar+r)) - sumx*log(xbar/(xbar+r));
                end
                
                % negative log-likelihood across all conditions
                nll = sum(nlls);
                
            end
        end
        
        function [nll,nlls] = negloglike_poisson(obj)
            % [nll,nlls] = negloglike_poisson(obj, tolBnd)
            %
            % Negative log-likihood for Poisson distribution fit

            nlls = zeros(obj.num_cond,1);
            
            for icond = 1:obj.num_cond
                % extract trials of this condition
                x = obj.spike_concat(1+obj.num_trial_cum(icond):obj.num_trial_cum(icond+1));
                
                % get sum, length and mean of x
                sumx = sum(x);
                n = length(x);
                xbar = mean(x);
                
                nlls(icond) = + n*xbar - sumx*log(xbar)+ sum(log(gamma(x+1)));
            end
            
            % negative log-likelihood across all conditions
            nll = sum(nlls);
            
       end
            
        function [pmf_table,fighandle] = pmf(obj,xval,showplot)
            % [pmf_table,fighandle] = pmf(obj,xval,showplot)
            % 
            % evaluate probability mass function based on negative binomial
            % parameter fit.
            %
            % Parameters
            % ----------
            % obj : object
            %   GainVariability class object
            % xval : array of int (optional)
            %   array of integers indicating at which x values the pmf
            %   should be evaluated. Default is xval=1:100;
            % showplot : bool (optional)
            %   boolean indicating whether to plot the pmf, default is
            %   false.
            %
            % Returns
            % -------
            % pmf_table : table
            %   table with y (nbinom) and y_poiss (Poisson) evaluated at x
            %   values.
            % fighandle : figure handle
            %
            
            % check input.
            if nargin<2 || isempty(xval)
                xval = 1:100;
            end
            if nargin<3
                showplot=false;
            end
            if isempty(obj.rhat)
                obj = obj.fit; % fit negative binomial
            end
            
            % compute pmf
            xbar = mean(obj.spike_concat); 
            r = 1/obj.rhat;
            p = r/(xbar+r);
            y = nbinpdf(xval,r,p); % fit negative binomial distribution
            y_poiss = poisspdf(xval,xbar); % fit poisson distribution
            
            % store in table
            pmf_table = table(xval(:),y(:),y_poiss(:),'VariableNames',{'x','y','y_poiss'});

            if showplot
                % plot both distributions on top of histogram
                figure(1),clf
                hold on
                histogram(obj.spike_concat,'Normalization','pdf')
                h(1) = plot(xval,y);
                h(2) = plot(xval,y_poiss);
                legend(h,{'nbinom','Poisson'});
                fighandle = gcf;
            else
                fighandle = [];
            end
        end
        
    end
end

