function selectivity = get_unit_selectivity(unitlist, type, subselection, threshold)


% perform checks on input
assert(istable(unitlist), 'unitlist should be a table')

if nargin<4
    threshold = 0.05;
end
if nargin<3
    subselection = [];
end

selectivity = false(height(unitlist),1);

% define selectivity - attention
cond1 = unitlist.s_att < threshold;
cond2 = unitlist.("s_att*dir") < threshold;
cond3 = unitlist.("s_att*dru") < threshold;
cond4 = unitlist.("s_att*dru*dir") < threshold;
selectivity_attention = cond1 | cond2 | cond3 | cond4;
clear cond*

% define selectivity - drug
cond1 = unitlist.s_dru < threshold;
cond2 = unitlist.("s_att*dru") < threshold;
cond3 = unitlist.("s_dru*dir") < threshold;
cond4 = unitlist.("s_att*dru*dir") < threshold;
selectivity_drug = cond1 | cond2 | cond3 | cond4;
clear cond*

% define selectivity - att*drug
cond1 = unitlist.("s_att*dru") < threshold;
cond2 = unitlist.("s_att*dru*dir") < threshold;
selectivity_att_drug = cond1 | cond2;
clear cond*

switch type
    case 'att'
        idx = selectivity_attention;
        
    case 'att&dru'
        idx = selectivity_attention & selectivity_drug;

    case 'att*dru'
        idx = selectivity_att_drug;
        
    case 'dru'
        idx = selectivity_drug;
        
    case 'visual'
        idx = unitlist.s_stim < threshold; % stim
        
    case 'visual&att'
        cond1 = unitlist.s_stim < threshold; % stim
        cond2 = selectivity_attention; % att

        idx = cond1 & cond2;

    case 'visual&att&dru'
        cond1 = unitlist.s_stim < threshold; % stim
        cond2 = selectivity_attention; % att
        cond3 = selectivity_drug; % dru
        
        idx = cond1 & cond2 & cond3;
        
    case 'visual&dru'
        cond1 = unitlist.s_stim < threshold; % stim
        cond2 = selectivity_drug; % dru

        idx = cond1 & cond2;
        
    case 'none'
        idx = true(height(P),1);
        
    otherwise
        error('unknown type')
end

selectivity(idx) = true;

%% subselection

if isempty(subselection)
    return
end

switch subselection{1}
    case 'subject'
        
    case 'drug'
        assert(iscell(subselection{2}), 'subselection{2} should be a cell')
        idx_subselection = unitlist.Drug==subselection{2};
        
end
assert(length(idx_subselection)==length(selectivity), 'subselection contains incorrect number of samples')

selectivity(~idx_subselection) = false;












