function [pstring,starstring] = get_significance_strings(p, varargin)
% convert p value into string for plotting
% [pstring,starstring] = get_significance_strings(p, varargin)
% 
% Parameters
% ----------
% p : array
%     array with p-value(s)
% varargin: cell, options include:
% 
%   - rounded: boolean, default 1, round p-value
%   - roundSet: '<' (default), 'power10'
%   - boldface: boolean, default 0, output string in bold face if p-value is < 0.05
%   - italic: boolean, default 0, output string in italic
%   - factorstring: string to indicate factor in analysis (e.g. p(attention))
%
% 
% Returns
% -------
% pstring : string 
%     string with p-value
% starstring : string 
%     string with stars indicating level of significance
%

get_varargin

if (~exist('rounded','var') || isempty(rounded)) && (~exist('round','var') || isempty(round))
    rounded = 1;
end

if ~exist('roundSet','var') || isempty(roundSet)
    roundSet = '<';
end

if ~exist('boldface','var') || isempty(boldface)
    boldface = 0;
end
if ~exist('italic','var') || isempty(italic)
    italic = 0;
end
if ~exist('factorstring','var') || isempty(factorstring)
    factorstring = [];
end


pstring = cell(length(p),1);
starstring = cell(length(p),1);
for ip = 1:length(p)
    if p(ip) >= 0.05
        if ~isempty(factorstring)
            pstring{ip} = sprintf('p(%s) = %.2f', factorstring{ip}, p(ip));
        else
            pstring{ip} = num2str(p(ip),'p = %.2f');
        end
        starstring{ip} = 'n.s.';
    elseif (p(ip) < 0.05) && (p(ip) >=0.01)
        if ~isempty(factorstring)
            pstring{ip} = sprintf('p(%s) < 0.05', factorstring{ip});
        else
            pstring{ip} = 'p < 0.05';
        end
        starstring{ip} = '*';
    elseif (p(ip) < 0.01) && (p(ip) >=0.001)
        if ~isempty(factorstring)
            pstring{ip} = sprintf('p(%s) < 0.01', factorstring{ip});
        else
            pstring{ip} = 'p < 0.01';
        end
        starstring{ip} = '**';
    elseif (p(ip) < 0.001)
        if ~isempty(factorstring)
            pstring{ip} = sprintf('p(%s) < 0.001', factorstring{ip});
        else
            pstring{ip} = 'p < 0.001';
        end
        starstring{ip} = '***';
    end
    
    %%% overwrite pstring in case of different rounding
    switch roundSet
        case 'power10'
            if (p(ip) == 0)
                pstring{ip} = sprintf('p = %d', p(ip));
            elseif (p(ip) >= 0.05)
                pstring{ip} = sprintf('p = %1.2f', p(ip));
            elseif (p(ip) < 0.05) && (p(ip) >= 0.01)
                %pstring{ip} = sprintf('p = %d', p(ip));
            elseif p(ip) < 0.01
                b = floor(log10(abs(p(ip)))) + 1;
                if ~isempty(factorstring)
                    pstring{ip} = sprintf('p(%s) < 10^{%d}', factorstring{ip}, b);
                else
                    pstring{ip} = sprintf('p < 10^{%d}', b);
                end
            end
    end
    
    if ~rounded
        
        if ~isempty(factorstring)
            pstring{ip} = sprintf('p(%s) = %.3f', factorstring{ip}, p(ip));
        else
            pstring{ip} = num2str(p(ip),'p = %.3f');
        end
        
        if (p(ip) < 0.001)
            if ~isempty(factorstring)
                pstring{ip} = sprintf('p(%s) < 0.001', factorstring{ip});
            else
                pstring{ip} = 'p < 0.001';
            end
        end
    end
       
    if boldface
        if p(ip) < 0.05
            pstring{ip} = ['\bf' pstring{ip}];
        end
    end
    
end

if length(pstring) == 1
    pstring = pstring{1};
    starstring = starstring{1};
end