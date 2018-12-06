function [ ret, IA, IB, IAB ] = wrap_dointersect( P, Q )
    global polytope_lib;
    
    if strcmp(polytope_lib, 'mpt') == 1
        [ret, IA, IB, IAB] = dointersect(P, Q);
    elseif strcmp(polytope_lib, 'pplmex') == 1
        %Qw = [];
        %Pw = [];
        %for i = 1:length(Q)
         %   Qw = [Qw polytope(Q(i).extreme)];
        %end
        %for i = 1:length(P)
         %   Pw = [Pw polytope(P(i).extreme)];
        %end
        %[ret, IA, IB, IAB] = dointersect(Pw, Qw);
        %return;
        intersect = zeros(length(P),length(Q));
        ret = 1;
        for i=1:length(P)
            for j=1:length(Q)
                if ~P(i).isdisjoint(Q(j))
                    intersect(i,j) = 1;
                else
                    ret = 0;
                end
            end
        end
        % TODO: IA, IB, IAB
    else
        error('Please define polytope_lib.');
    end

end

