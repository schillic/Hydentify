function [ ret ] = wrap_polytope( H, K )
    %WRAP_POLYTOPE Summary of this function goes here
    %   Detailed explanation goes here
    global polytope_lib;

    if strcmp(polytope_lib, 'mpt') == 1
        if nargin == 1
            ret = polytope(H);
        else
            ret = polytope(H, K);
        end
    elseif strcmp(polytope_lib, 'pplmex') == 1
        if nargin == 1
            ret = ppl_polytope(H);
        else
            ret = ppl_polytope(-H,K);
        end
    else
        error('Please define polytope_lib.');
    end
end

