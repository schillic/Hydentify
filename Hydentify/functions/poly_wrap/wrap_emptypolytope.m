function [ ret ] = wrap_emptypolytope()
%WRAP_EMPTYPOLYTOPE Summary of this function goes here
%   Detailed explanation goes here
    global polytope_lib;
    
    if strcmp(polytope_lib, 'mpt') == 1
        ret = polytope(); % the empty polytope
    elseif strcmp(polytope_lib, 'pplmex') == 1
        ret = ppl_polytope(0, -1);
    else
        error('Please define polytope_lib.');
    end

end

