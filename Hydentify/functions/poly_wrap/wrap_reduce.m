function [ ret ] = wrap_reduce( P )
%WRAP_REDUCE Summary of this function goes here
%   Detailed explanation goes here
    global polytope_lib;
    
    if strcmp(polytope_lib, 'mpt') == 1
        ret= reduce(P);
    elseif strcmp(polytope_lib, 'pplmex') == 1
        ret = P; % any 'real' operation should minimize the resulting polytope (feature of the PPL wrapper).s
    else
        error('Please define polytope_lib.');
    end

end

