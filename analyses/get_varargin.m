%%% unpack varargin into separate variables
%
% Jochem van Kempen. 19-06-2018


% get specific varargin
for iarg = 1:2:length(varargin)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% other
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmpi(varargin{iarg}, 'NA_action')
        NA_action = varargin{iarg+1};
        
    elseif strcmpi(varargin{iarg}, 'rounded') || strcmpi(varargin{iarg}, 'round') 
        rounded = varargin{iarg+1};
        
    elseif strcmpi(varargin{iarg}, 'roundSet')
        roundSet = varargin{iarg+1};
        
    elseif strcmpi(varargin{iarg}, 'factorstring')
        factorstring = varargin{iarg+1};
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% otherwise
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        warning('Unknown parameter %s, varargin', varargin{iarg})
        if isnumeric(varargin{iarg+1}) % numeric
            eval([varargin{iarg} ' = [' num2str(varargin{iarg+1}) '];'])
        elseif ischar(varargin{iarg+1}) % strings
            eval([varargin{iarg} ' = ''' varargin{iarg+1} ''';'])
        elseif iscell(varargin{iarg+1}) % cell arrays, with strings
            nel = length(varargin{iarg+1});
            eval([varargin{iarg} ' = cell(1, ' num2str(nel) ') ;'])
            for iel = 1:nel
                eval([varargin{iarg} '{iel} = ''' varargin{iarg+1}{iel} ''';'])
            end
        else
            error('Cannot extract parameter, varargin')
        end
        
    end
    
    
    
end

