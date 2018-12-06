function [ V, R, Pb ] = wrap_extreme( P )
%WRAP_EXTREME Summary of this function goes here
%   Detailed explanation goes here
    global polytope_lib;
    
    if strcmp(polytope_lib, 'mpt') == 1
        [V,R,Pb] = extreme(P);
    elseif strcmp(polytope_lib, 'pplmex') == 1
        V = P.extreme;
        R = [];
        Pb = P;
    else
        error('Please define polytope_lib.');
    end

end

