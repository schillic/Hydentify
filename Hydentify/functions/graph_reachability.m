function [reached] = graph_reachability(transition_relation, initial_boxes, num_boxes, print_result, reachable_locations)
%GRAPH_REACHABILITY Computes the reachable boxes in the SpaceEx model
% formula for the initial box:
% Box_x1_x2_..._xn translates to the index
% I = x1 + sum_i>1 [(xi - 1) prod_j=1^i-1 |vj|]
% where |vj| is the number of intervals in dimension j
%
% -------------------------------------------------------------------------
% Author:  Christian Schilling
%          University of Freiburg, Germany
% email:   schillic at informatik.uni-freiburg.de
% Website: http://swt.informatik.uni-freiburg.de/staff/christian_schilling
%
% September 2015
% -------------------------------------------------------------------------

    global g_benchmark_reachable_locations_graph_reachability;
    
    % set up the data structure
    src2dest = cell(num_boxes, 1);
    src = transition_relation(:, 1);
    dest = transition_relation(:, 2);
    for i = 1 : size(src)
       if (src(i) ~= dest(i)) % skip the self-loop
           src2dest{src(i)}(length(src2dest{src(i)}) + 1) = dest(i);
       end
    end
    reached = zeros(length(src2dest), 1);
    
    % initialize stack with the initial boxes
    stack = zeros(length(src2dest), 1);
    stack_size = size(initial_boxes, 2);
    for i = 1 : stack_size
        stack(i) = initial_boxes(1, i);
        reached(stack(i)) = 1;
    end
    
    if (~ isempty(reachable_locations))
        reached = markReachableLocations(reachable_locations, reached);
    end
    
    % depth-first search
    while (stack_size > 0)
        src = stack(stack_size);
        stack_size = stack_size - 1;

        % visit successors
        targets = src2dest{src};
        for t = targets
            if (reached(t) == 0)
                stack_size = stack_size + 1;
                stack(stack_size) = t;
                reached(t) = 1;
            end
        end
    end
    
    if (print_result)
        count = 0;
        for i = 1 : length(reached)
            if (reached(i) == 1)
                count = count + 1;
            end
        end
        fprintf('\ngraph reachability ');
        if (~ isempty(reachable_locations))
            fprintf('(+ previous analysis) ');
        end
        fprintf('saves %d/%d locations\n', num_boxes - count, num_boxes);
        g_benchmark_reachable_locations_graph_reachability = count;
    end
end


function [ reached ] = markReachableLocations(reachable_locations, reached)
% marks the reachable locations from a previous iteration
% has a 1-entry iff the location is not reachable

    i_next_reachable = 0;
    do_continue = 1;
    for i = 1 : length(reached)
        if (do_continue == 1)
            if (i_next_reachable == length(reachable_locations))
                i_next = length(reached) + 1;
                do_continue = -1;
            else
                i_next_reachable = i_next_reachable + 1;
                i_next = reachable_locations{i_next_reachable};
            end
        end
        
        if (i < i_next)
            % do not remove initial locations
            % Reason: It would cause additional efforts to filter them from the
            % configuration.
            if (reached(i) == 0)
                reached(i) = -1;
            end
            do_continue = 0;
        else
            assert(i == i_next, ...
                'The entries should be sorted in ascending order.');
            if (do_continue == 0)
                do_continue = 1;
            end
        end
    end
end