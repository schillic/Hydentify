function [ R, kept ] = wrap_reduceunion( P )
%WRAP_REDUCEUNION Summary of this function goes here
%   Detailed explanation goes here
    global polytope_lib;
    
    if strcmp(polytope_lib, 'mpt') == 1
        [R, kept] = reduceunion(P);
    elseif strcmp(polytope_lib, 'pplmex') == 1
        kept = ones(1,length(P));
        
        for i=1:length(P)
            for j=1:length(P)
                if i~=j
                    if kept(i) && P(i).contains(P(j))
                        kept(j) = 0; % i contains j, don't need i
                    end
                end
            end
        end
        
        %for i=1:length(P)
        %    for j=1:length(P)
        %       if i~=j
        %           tmp = P(i) | P(j);
        %            if ~(kept(i) && tmp ~= P(i))
        %                kept(i) = 0;
        %            end
        %       end
        %   end
        %end
        R = [];
        for i=1:length(P)
            if kept(i)
                R = [R, P(i)];
            end
        end
    else
        error('Please define polytope_lib.');
    end

end

