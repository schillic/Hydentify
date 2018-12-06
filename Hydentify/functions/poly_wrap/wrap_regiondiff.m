function [ R ] = wrap_regiondiff( P, Pb, options )
%WRAP_REDUCEUNION Summary of this function goes here
%   Detailed explanation goes here
    global polytope_lib;
    
    if strcmp(polytope_lib, 'mpt') == 1
        R=regiondiff(P, Pb, options);
    elseif strcmp(polytope_lib, 'pplmex') == 1
        Pbt = [];
        R = [];
        for i=1:length(Pb)
            Pbt = [Pbt polytope(Pb(i).extreme)];
        end
        tmp=regiondiff(polytope(P.extreme), Pbt, options);
        for i = 1:length(tmp)
            if ~isfulldim(tmp(i))
                R = [R ppl_polytope(0, -1)];
            else
                R = [R ppl_polytope(extreme(tmp(i)))];
            end
        end
    end
end

