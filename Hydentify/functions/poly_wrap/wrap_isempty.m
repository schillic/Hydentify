function [ ret ] = wrap_isempty( P )
%WRAP_ISEMPTY Summary of this function goes here
%   Detailed explanation goes here
    global polytope_lib;
    
    if strcmp(polytope_lib, 'mpt') == 1
        ret= ~isfulldim(P);
    elseif strcmp(polytope_lib, 'pplmex') == 1
        if length(P) > 1
            ret = zeros(length(P),1);
            for i=1:length(P)
                if(P(i).is_empty)
                    ret(i) = 1;
                end
            end
        elseif length(P) == 0
            ret = 1;
        else
            ret = P.is_empty;
        end
    else
        error('Please define polytope_lib.');
    end

end

