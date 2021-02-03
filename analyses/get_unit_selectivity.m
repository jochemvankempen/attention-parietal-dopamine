function selectivity = get_unit_selectivity(P, type, threshold)

assert(istable(P), 'P is not a table')

if nargin<3
    threshold = 0.05;
end

selectivity = false(height(P),1);

switch type
    
    case 'att'
        
        cond1 = P.att < threshold;
        cond2 = P.("att*dir") < threshold;
        cond3 = P.("att*dru") < threshold;
        cond4 = P.("att*dru*dir") < threshold;
        
        idx = cond1 | cond2 | cond3 | cond4;
        
    case 'att+dru'

        cond1 = P.("att") < threshold;
        cond2 = P.("dru") < threshold;

        idx = cond1 & cond2;

    case 'att*dru'

        cond1 = P.("att*dru") < threshold;

        idx = cond1;
end

selectivity(idx) = true;

