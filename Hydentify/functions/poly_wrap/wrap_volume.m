function V = wrap_volume( P, options )
%WRAP_VOLUME Summary of this function goes here
%   Detailed explanation goes here
    global polytope_lib;
    
    if strcmp(polytope_lib, 'mpt') == 1
        V=volume(P);
    elseif strcmp(polytope_lib, 'pplmex') == 1
        P1 = polytope(P.extreme);
        V = volume(P1);
    end
end

