function analyze_parametric_transition_system_hydentify_iterative()
%ANALYZE_PARAMETRIC_TRANSITION_SYSTEM_HYDENTIFY_ITERATIVE Iterative analysis in Hydentify.
% This function allows for more specified analysis of nodes in the search tree
% in several fashions.
% TODO detailed explanation
%
% -------------------------------------------------------------------------
%
% Many parts taken from 'analyze_parametric_transition_system_hydentify.m'.
%
% Authors: Christian Schilling
%          University of Freiburg, Germany
% email:   schillic at informatik.uni-freiburg.de
% Website: http://swt.informatik.uni-freiburg.de/staff/christian_schilling
%
% September 2015
% -------------------------------------------------------------------------

global analysis_type;
global Pspace;
global constraint_list;
global live_state;
global result;
global transition_relation_existsP;
global transition_relation_forallP;
global transient_set_existsP;
global transient_set_forallP;
global G_DEBUG_OUTPUT;
global G_PRINT_INTERMEDIATE_RESULT;
global G_ANALYSIS_STRATEGY;
global G_SPACEEX_HEURISTICS;
global G_NODE_HEURISTICS;
global G_STORE_SEARCH_TREE;
global G_HYPY_BOUND;

% benchmark data
global g_benchmark_rovergene_exists_count;
global g_benchmark_rovergene_forall_count;
global g_benchmark_hydentify_exists_count;
global g_benchmark_hydentify_forall_count;
global g_benchmark_rovergene_exists_time;
global g_benchmark_rovergene_forall_time;
global g_benchmark_hydentify_exists_time;
global g_benchmark_hydentify_forall_time;
global g_benchmark_total_time;
global g_benchmark_total_parameter_space_size;
global g_benchmark_verified_parameter_space_size;
global g_benchmark_verified_parameter_space_percentage;
global g_benchmark_nodes_explored;
global g_benchmark_reachable_locations_graph_reachability;
global g_benchmark_reachable_locations_spaceex_reachability_total;
g_benchmark_rovergene_exists_count = 0;
g_benchmark_rovergene_forall_count = 0;
g_benchmark_hydentify_exists_count = 0;
g_benchmark_hydentify_forall_count = 0;
g_benchmark_rovergene_exists_time = 0;
g_benchmark_rovergene_forall_time = 0;
g_benchmark_hydentify_exists_time = 0;
g_benchmark_hydentify_forall_time = 0;
g_benchmark_total_parameter_space_size = 0;
g_benchmark_verified_parameter_space_size = 0;
g_benchmark_verified_parameter_space_percentage = 0;
g_benchmark_nodes_explored = 0;
g_benchmark_total_time = 0;
g_benchmark_reachable_locations_graph_reachability = 0;
g_benchmark_reachable_locations_spaceex_reachability_total = 0;

% number of constraints
CONSTRAINTS_LENGTH = length(constraint_list.K);

if (G_DEBUG_OUTPUT)
    % print constraints
    fprintf('constraints:\n');
    for i = 1 : CONSTRAINTS_LENGTH
        row = constraint_list.H(i, :);
        for j = 1 : size(row, 2)
            if (j > 1)
                fprintf(' + %d p%d', row(j), j);
            else
                fprintf('   %d p%d', row(j), j);
            end
        end
        fprintf(' <= %d\n', constraint_list.K(i));
    end
end

if (G_STORE_SEARCH_TREE)
    if (CONSTRAINTS_LENGTH < 7)
        INITIAL_TREE_SIZE = 2^(CONSTRAINTS_LENGTH + 1) - 1;
    else
        INITIAL_TREE_SIZE = 100;
    end
end

% true iff parameter space is empty
is_Pspace_empty = wrap_isempty(Pspace);
if (is_Pspace_empty)
    % Create a dummy stack with one frame.
    stack = cell(1, 1);
else
    % precompute parameter space representation
    [PbH_PURE, PbK_PURE]= wrap_hk(Pspace);
    [POLYTOPE_HEIGHT, POLYTOPE_WIDTH] = size(PbH_PURE);
    
    % Create a stack for depth-first traversal of the search tree.
    % It has size at most (#constraints + 1).
    stack = cell(CONSTRAINTS_LENGTH + 1, 1);
end

% initialize with root node
stack_size = 1;
stack{1} = [];
switch (G_ANALYSIS_STRATEGY)
    case {EOptions.ANALYSIS_ROVERGENE, EOptions.ANALYSIS_SPACEEX_ONLY}
        % nothing to do
        
    case {EOptions.ANALYSIS_EAGER, EOptions.ANALYSIS_LAZY}
        % additional stack for SpaceEx calls
        stack_size_spaceex = 0;
        stack_spaceex = cell(length(stack), 2);

    otherwise
        assert(false, 'Unknown strategy.');
end

% HyPy: additional map for reachable locations
if (G_HYPY_BOUND >= 0)
    map_hypy = containers.Map({'-1', '[]'}, {{}, {}});
else
    reachable_locations = [];
    spaceex_mode = EOptions.SPACEEX_EARLY_TERMINATION;
end

if (G_STORE_SEARCH_TREE)
    % explicit search tree in memory
    tree = cell(INITIAL_TREE_SIZE, 7);
    tree_size = 0;
    add_tree_node(1, 0, []);
else
    % dummy value
    node_idx = 0;
end

% Boolean for remembering the consistency test result of the last lhs sibling.
% This way redundancy tests only have to be performed for the rhs siblings.
is_last_inconsistency_test_result = false;

% traverse search tree
while (true)
    % check termination condition (depends on analysis mode)
    switch (G_ANALYSIS_STRATEGY)
        case {EOptions.ANALYSIS_ROVERGENE, EOptions.ANALYSIS_SPACEEX_ONLY}
            % one stack
            if (stack_size == 0)
                break;
            end
        
        case {EOptions.ANALYSIS_EAGER, EOptions.ANALYSIS_LAZY}
            % two stacks
            if (stack_size + stack_size_spaceex == 0)
                break;
            end
        
        otherwise
            assert(false, 'Unknown strategy.');
    end
    
    % decide on which search node to continue with
    switch (G_ANALYSIS_STRATEGY)
        case {EOptions.ANALYSIS_ROVERGENE, EOptions.ANALYSIS_SPACEEX_ONLY}
            % choose from the only stack
            b = pop_from_normal_stack();
            CURRENT_MODE = EOptions.MODE_NORMAL;
        
        case EOptions.ANALYSIS_EAGER
            % prefer SpaceEx (does not matter, but is pure depth-first search)
            if (stack_size_spaceex == 0)
                b = pop_from_normal_stack();
                CURRENT_MODE = EOptions.MODE_NORMAL;
            else
                [b, CURRENT_MODE] = pop_from_spaceex_stack();
            end
        
        case EOptions.ANALYSIS_LAZY
            % choose from two possible stacks
            switch (G_SPACEEX_HEURISTICS)
                case EOptions.LAZY_SPACEEX_FIRST
                    if (stack_size_spaceex == 0)
                        b = pop_from_normal_stack();
                        CURRENT_MODE = EOptions.MODE_NORMAL;
                    else
                        [b, CURRENT_MODE] = pop_from_spaceex_stack();
                        % TODO With this strategy ("prefer SpaceEx") we should
                        % probably search through the regular stack in a
                        % successful exists case to avoid RoVerGeNe analysis of
                        % the other child node (if any). A constant lookup
                        % should be enough (it would be the topmost element).
                    end
                
                case EOptions.LAZY_SPACEEX_LAST
                    if (stack_size == 0)
                        [b, CURRENT_MODE] = pop_from_spaceex_stack();
                    else
                        b = pop_from_normal_stack();
                        CURRENT_MODE = EOptions.MODE_NORMAL;
                    end
            end
        
        otherwise
            assert(false, 'Unknown strategy.');
    end
    
    % 'b' is a Boolean vector encoding the constraints
    B_LENGTH = length(b);
    
    if (G_DEBUG_OUTPUT)
        fprintf('\n-------------------\n');
        if (isempty(b))
            fprintf('\nb = [] \n');
        else
            fprintf('\nb ='), disp(b);
        end
    end
    
    % use HyPy for the first n levels in the search hierarchy
    if (G_HYPY_BOUND >= 0)
        if (B_LENGTH <= G_HYPY_BOUND)
            spaceex_mode = EOptions.SPACEEX_HYPY;
            
            % reachable locations (HyPy feature; empty = all locations)
            bit_string = get_bit_string(b);
            assert(map_hypy.isKey(bit_string), 'The key must be contained.');
            reachable_locations = map_hypy(bit_string);
            map_hypy.remove(bit_string);
        else
            spaceex_mode = EOptions.SPACEEX_EARLY_TERMINATION;
        end
    end
    
    if (G_STORE_SEARCH_TREE)
        % find search tree node
        node_idx = 1;
        for idx = b
            if (idx)
                node_idx = tree{node_idx, EOptions.TREE_IDX_RIGHT};
            else
                node_idx = tree{node_idx, EOptions.TREE_IDX_LEFT};
            end
        end
        
        % reset stack entry
        tree{node_idx}(EOptions.TREE_IDX_STACK) = 0;
    end
    
    % This Boolean is set to 1 if we detect that the newly added constraint
    % is inconsistent with the previous ones.
    is_constraint_inconsistent = false;
    
    % This Boolean is set to 1 if we detect that the newly added constraint
    % is redundant with the previous ones. In that case, the two
    % corresponding parameter sets are equal and we know that the previous
    % tests yielded answer_existsP = 0 and answer_forallP = 1 (the only case
    % where we continue the computations). So we can set directly the same
    % values and bypass the verification tests.
    is_constraint_redundant = false;
    
    % create the polytope from vector
    if (is_Pspace_empty)
        Pb = wrap_emptypolytope();
    elseif (B_LENGTH == 0)
        PbH = cat(1, PbH_PURE);
        PbK = cat(1, PbK_PURE);
        Pb = wrap_polytope(PbH, PbK);
        
        if (G_STORE_SEARCH_TREE)
            tree{node_idx, EOptions.TREE_IDX_SKIP} = EOptions.SKIP_NORMAL;
        end
        
        if (G_DEBUG_OUTPUT)
            fprintf('root node, skipping redundancy and inconsistency tests\n');
        end
    else
        % optimized: no reallocation all the time
        PbH = cat(1, PbH_PURE, zeros(B_LENGTH, POLYTOPE_WIDTH));
        PbK = cat(1, PbK_PURE, zeros(B_LENGTH, 1));
        j = POLYTOPE_HEIGHT;
        for i = 1 : B_LENGTH
            j = j + 1;
            if b(i)
                PbH(j, :) = constraint_list.H(i, :);
                PbK(j) = constraint_list.K(i);
            else
                PbH(j, :) = -constraint_list.H(i, :);
                PbK(j) = -constraint_list.K(i);
            end
        end
        Pb = wrap_polytope(PbH, PbK);
        
        % skip when adding constraints made parameter space empty
        if (b(B_LENGTH))
            % an rhs sibling on the stack is never inconsistent
            is_skip_emptiness_test = true;
        else
            is_skip_emptiness_test = false;
            switch (CURRENT_MODE)
                case {EOptions.MODE_SPACEEX_EX}
                    % nothing to do, as we know the parent was not inconsistent
                    is_skip_emptiness_test = true;
                
                case {EOptions.MODE_NORMAL, EOptions.MODE_SPACEEX_BOTH, ...
                      EOptions.MODE_SPACEEX_BOTH_FIRST}
                    if (G_STORE_SEARCH_TREE)
                        % The result can be derived from the tree iff the node
                        % or the sibling was already visited.
                        switch (tree{node_idx, EOptions.TREE_IDX_SKIP})
                            case {EOptions.SKIP_INCONSISTENT, ...
                                  EOptions.SKIP_INCONSISTENT_INFERRED}
                                % already known, inconsistent
                                is_constraint_inconsistent = true;
                                is_skip_emptiness_test = true;

                            case 0
                                assert(~ b(B_LENGTH), ...
                                    'The node should be an lhs child.');

                            case {EOptions.SKIP_NORMAL, ...
                                  EOptions.SKIP_NORMAL_INFERRED, ...
                                  EOptions.SKIP_REDUNDANT, ...
                                  EOptions.SKIP_REDUNDANT_INFERRED}
                                % already known, not inconsistent
                                is_skip_emptiness_test = true;

                            otherwise
                                assert(false, 'Unknown mode.');
                        end
                    end

                otherwise
                    assert(false, 'Unknown mode.');
            end
        end
        
        if (~ is_skip_emptiness_test)
            % perform emptiness test on polytope
            is_constraint_inconsistent = wrap_isempty(Pb);
            
            if (is_constraint_inconsistent)
                if (G_STORE_SEARCH_TREE)
                    % update tree data
                    tree{node_idx, EOptions.TREE_IDX_SKIP} = ...
                        EOptions.SKIP_INCONSISTENT;
                    if (~ b(B_LENGTH))
                        % infer redundancy status for rhs sibling
                        tree{node_idx + 1, EOptions.TREE_IDX_SKIP} = ...
                            EOptions.SKIP_REDUNDANT_INFERRED;
                    end
                end
            end
        else
            if (G_DEBUG_OUTPUT)
                fprintf('skipping emptiness test\n');
            end
        end
        
        if (is_constraint_inconsistent)
            assert(B_LENGTH > 0, 'The root node cannot be inconsistent.');
            
            % set flag for skipping redundancy test for rhs sibling
            if (~ b(B_LENGTH))
                is_last_inconsistency_test_result = true;
            end
            
            % skip current node
            if (G_DEBUG_OUTPUT)
                fprintf('current constraints are inconsistent, skipping\n');
            end
            
            % possibly push the transitive parent (see function description)
            switch (CURRENT_MODE)
                case EOptions.MODE_NORMAL
                    check_and_push_remember_node();

                case {EOptions.MODE_SPACEEX_EX, ...
                      EOptions.MODE_SPACEEX_BOTH, ...
                      EOptions.MODE_SPACEEX_BOTH_FIRST}
                    % nothing to do

                otherwise
                    assert(false, 'Unknown mode.');
            end
            
            continue;
        end
        
        if (b(B_LENGTH))
            % for rhs siblings the redundancy is inferred
            is_constraint_redundant = is_last_inconsistency_test_result;
            is_last_inconsistency_test_result = false;
            is_skip_redundancy_test = true;
        else
            is_skip_redundancy_test = false;
            if (G_STORE_SEARCH_TREE)
                % The result can be derived from the tree iff the node or the
                % sibling was already visited.
                switch (tree{node_idx, EOptions.TREE_IDX_SKIP})
                    case {EOptions.SKIP_REDUNDANT, ...
                          EOptions.SKIP_REDUNDANT_INFERRED}
                        % already known, redundant
                        is_constraint_redundant = true;
                        is_skip_redundancy_test = true;

                    case 0
                        if (B_LENGTH > 0)
                            assert(~ b(B_LENGTH), ...
                                'The node should be an lhs child.');
                        end

                    case {EOptions.SKIP_NORMAL, EOptions.SKIP_NORMAL_INFERRED}
                        % already known, not inconsistent
                        is_skip_redundancy_test = true;

                    case {EOptions.SKIP_INCONSISTENT, ...
                          EOptions.SKIP_INCONSISTENT_INFERRED}
                        % already known, not inconsistent
                        is_skip_redundancy_test = true;

                        % inconsistent nodes should have been skipped before
                        assert(false, 'unreachable code');

                    otherwise
                        assert(false, 'Unknown mode.');
                end
            end
        end
        
        if (~ is_skip_redundancy_test)
            % We test redundancy of newly added constraint wrt. the previous
            % ones.
            % We use linear programming as described in FAQ on polyhedral
            % computation (Fukuda, 2004), sec. 2.21.
            % We have PbH(last), PbK(last) is redundant iff f* <= PbK(last)
            % with f* = max PbH(last) p, s.t. PbH p <= B, where B = PbK
            % except that B(last) = PbK(last) + 1.
            % Also, we use [xopt, fval] = mpt_solveLP(-PbK(last), PbH, B).
            % Note that we use -h because mpt_solve minimizes.
            B = PbK;
            B(length(PbK)) = PbK(length(PbK)) + 1;
            [~, fval] = mpt_solveLP(-PbH(length(PbK), :), PbH, B);
            if (-fval <= PbK(length(PbK)))
                % the added constraint is redundant
                % (-fval because we used minimization)
                is_constraint_redundant = true;
            end
            
            if (G_STORE_SEARCH_TREE)
                if (is_constraint_redundant)
                    % update tree data
                    tree{node_idx, EOptions.TREE_IDX_SKIP} = ...
                        EOptions.SKIP_REDUNDANT;
                else
                    tree{node_idx, EOptions.TREE_IDX_SKIP} = ...
                        EOptions.SKIP_NORMAL;
                    if (~ b(B_LENGTH))
                        % infer normal status for rhs sibling
                        tree{node_idx + 1, EOptions.TREE_IDX_SKIP} = ...
                            EOptions.SKIP_NORMAL_INFERRED;
                    end
                end
            end
        else
            if (G_DEBUG_OUTPUT)
                fprintf('skipping redundancy test\n');
            end
        end
    end
    
    % true iff search tree should be descended (refined later)
    % necessary condition: analysis type is parameter synthesis
    
    is_descend_in_search_tree = ((analysis_type == 1) && ...
                                (B_LENGTH < CONSTRAINTS_LENGTH));
    
    if (is_constraint_redundant)
        assert(B_LENGTH > 0, 'The root node cannot be redundant.');
        
        if (G_DEBUG_OUTPUT)
            fprintf('new constraint is redundant, skipping\n');
        end
        
        switch (CURRENT_MODE)
            case EOptions.MODE_NORMAL
                % nothing to do, just update live_states
                live_state.existsP{B_LENGTH + 1} = live_state.existsP{B_LENGTH};
                live_state.forallP{B_LENGTH + 1} = live_state.forallP{B_LENGTH};
                
                % optimization: remove sibling from stack
                % The idea is that for two siblings s1, s2 the following holds:
                %   s1 is inconsistent iff s2 is redundant
                % So when we find a redundant node, we can just remove the
                % sibling as it will be inconsistent and so we would skip it
                % anyway.
                if (~ b(B_LENGTH))
                    % This is an lhs sibling, i.e.,
                    % the rhs sibling is the topmost element on the stack.
                    % As we know it is inconsistent, we just pop it.
                    stack_size = stack_size - 1;
                    
                    if (G_STORE_SEARCH_TREE)
                        assert(...
                            tree{node_idx + 1, EOptions.TREE_IDX_STACK} == ...
                                stack_size + 1, ...
                            'The node must be the former topmost stack element.');
                        tree{node_idx + 1, EOptions.TREE_IDX_STACK} = 0;
                        tree{node_idx + 1, EOptions.TREE_IDX_SKIP} = ...
                            EOptions.SKIP_INCONSISTENT_INFERRED;
                    end

                    if (G_DEBUG_OUTPUT)
                        fprintf('sibling constraints are inconsistent, popping\n');
                    end
                end
                
                % possibly push the transitive parent (see function description)
                check_and_push_remember_node();
            
            case EOptions.MODE_SPACEEX_EX
                % do not descend but ascend
                is_descend_in_search_tree = false;
                push_to_spaceex_stack(b, true, CURRENT_MODE, true, node_idx);
            
            case {EOptions.MODE_SPACEEX_BOTH, EOptions.MODE_SPACEEX_BOTH_FIRST}
                % optimization: remove rhs sibling from stack (s. above)
                if (~ b(B_LENGTH))
                    stack_size_spaceex = stack_size_spaceex - 1;
                    
                    if (G_STORE_SEARCH_TREE)
                        tree{node_idx, EOptions.TREE_IDX_SKIP} = ...
                            EOptions.SKIP_REDUNDANT;
                        tree{node_idx + 1, EOptions.TREE_IDX_SKIP} = ...
                            EOptions.SKIP_INCONSISTENT_INFERRED;
                    end

                    if (G_DEBUG_OUTPUT)
                        fprintf('sibling constraints are inconsistent, popping\n');
                    end
                end
            
            otherwise
                    assert(false, 'Unknown mode.');
        end
    else
        % extreme points
        [~, ~, Pb] = wrap_extreme(Pb);
        
        if (G_DEBUG_OUTPUT)
            fprintf('Computing concrete transition system: ... ');
        end
        
        % This variable holds the last node explored by RoVerGeNe.
        % If in the lazy strategy all children are redundant/inconsistent, the
        % current node should be analyzed by SpaceEx. This is forgotten here, as
        % we have found another child which can be explored.
        remember_node = -1;
        remember_node_idx = -1;
        
        % compute and store the transition relation for existential and
        % universal semantics
        [transition_relation_existsP, transition_relation_forallP] = ...
            compute_transition_relation(Pb);
        
        if (G_DEBUG_OUTPUT)
            fprintf('|%g s|\n', toc);
        end
        
        switch (CURRENT_MODE)
            case EOptions.MODE_NORMAL
                if (G_DEBUG_OUTPUT)
                    fprintf('Computing transient sets: ... ');
                end
                
                % compute and store the set of transient states for existential
                % and universal semantics (only necessary for RoVerGeNe)
                [transient_set_existsP, transient_set_forallP] = ...
                    compute_transient_set(b, Pb, transition_relation_existsP, ...
                                            transition_relation_forallP);
                
                if (G_DEBUG_OUTPUT)
                    fprintf('|%g s|\n', toc);
                end
            
            case {EOptions.MODE_SPACEEX_EX, EOptions.MODE_SPACEEX_BOTH, ...
                  EOptions.MODE_SPACEEX_BOTH_FIRST}
                % nothing to do
            
            otherwise
                assert(false, 'Unknown mode.');
        end
        
        % ---------------
        % actual analysis
        %
        % Depending on the global option, there are several strategies to use:
        %
        % 1) eager approach (see strategy_eager()):
        % Use SpaceEx all the time RoVerGeNe was not good enough (the
        % "classical" approach).
        %
        % 2) lazy approach (see strategy_lazy()):
        % Use SpaceEx only after RoVerGeNe found a valid parameter set or wants
        % to prune the search tree.
        % ---------------
        switch (G_ANALYSIS_STRATEGY)
            case EOptions.ANALYSIS_ROVERGENE
                [is_valid, is_descend_in_search_tree] = ...
                    strategy_rovergene(is_descend_in_search_tree);
            
            case EOptions.ANALYSIS_SPACEEX_ONLY
                [is_valid, is_descend_in_search_tree, reachable_locations] = ...
                    strategy_spaceex_only(is_descend_in_search_tree, ...
                        reachable_locations, spaceex_mode);
            
            case EOptions.ANALYSIS_EAGER
                [is_valid, is_descend_in_search_tree, reachable_locations] = ...
                    strategy_eager(is_descend_in_search_tree, ...
                    reachable_locations, spaceex_mode);
            
            case EOptions.ANALYSIS_LAZY
                [is_valid, is_descend_in_search_tree, reachable_locations] = ...
                    strategy_lazy(is_descend_in_search_tree, ...
                    reachable_locations, spaceex_mode);
            
            otherwise
                assert(false, 'Unknown strategy.');
        end
        
        % store a valid parameter set
        if (is_valid)
            result.property_robustly_satisfied = 2;
            result.valid_parameter_set.certain = ...
                [result.valid_parameter_set.certain, Pb];
            
            if (G_PRINT_INTERMEDIATE_RESULT)
                fprintf('intermediate result: ');
                display_result();
            end
        end
    end
    
    % assert tree information validity
    if (G_STORE_SEARCH_TREE)
        switch (tree{node_idx, EOptions.TREE_IDX_SKIP})
            case {EOptions.SKIP_NORMAL, EOptions.SKIP_NORMAL_INFERRED}
                switch (CURRENT_MODE)
                    case EOptions.MODE_NORMAL
                        assert((tree{node_idx, EOptions.TREE_IDX_DATA}...
                            (EOptions.DATA_IDX_RGE) ~= 0), ...
                            'there should be some result for RG exists analysis');

                    case {EOptions.MODE_SPACEEX_EX, ...
                          EOptions.MODE_SPACEEX_BOTH, ...
                          EOptions.MODE_SPACEEX_BOTH_FIRST}
                        assert((tree{node_idx, EOptions.TREE_IDX_DATA}...
                                (EOptions.DATA_IDX_SRE) ~= 0), ...
                            'there should be some result for RG exists analysis');

                    otherwise
                        assert(false, 'Unknown mode.');
                end
            
            case {EOptions.SKIP_REDUNDANT, EOptions.SKIP_REDUNDANT_INFERRED, ...
                  EOptions.SKIP_INCONSISTENT, ...
                  EOptions.SKIP_INCONSISTENT_INFERRED}
                % do nothing
            
            otherwise
                assert(false, 'Unknown mode.');
        end
    end
    
    % push new nodes on the stack (depending on the mode)
    if (is_descend_in_search_tree)
        switch (CURRENT_MODE)
            case EOptions.MODE_NORMAL
                % push two child nodes to the normal stack to continue search on
                % the next level
                for next_bit = 1 : -1 : 0
                    push_to_normal_stack([b next_bit]);
                end
                
                if (G_STORE_SEARCH_TREE)
                    % add nodes to search tree in correct order
                    add_tree_node(stack_size, node_idx, [b 0]);
                    add_tree_node(stack_size - 1, -node_idx, [b 1]);
                end
            
            case EOptions.MODE_SPACEEX_EX
                assert(false, 'unreachable code');
            
            case {EOptions.MODE_SPACEEX_BOTH, EOptions.MODE_SPACEEX_BOTH_FIRST}
                % push two child nodes to the SpaceEx analysis stack to continue
                % search on the next level
                for next_bit = 1 : -1 : 0
                    push_to_spaceex_stack([b next_bit], false, ...
                                  EOptions.MODE_SPACEEX_BOTH, false, node_idx);
                end
                
                if (G_STORE_SEARCH_TREE)
                    % add nodes to search tree in correct order
                    add_tree_node(-stack_size_spaceex, node_idx, [b 0]);
                    add_tree_node(-stack_size_spaceex + 1, -node_idx, [b 1]);
                end
            
            otherwise
                assert(false, 'Unknown mode.');
        end
        
        % update reachable locations
        if (G_HYPY_BOUND > 0)
            for next_bit = 1 : -1 : 0
                map_hypy(get_bit_string([b next_bit])) = reachable_locations;
            end
        end
    end
end

% - write final benchmark data -

% about valid parameter set size
[g_benchmark_total_parameter_space_size, ...
    g_benchmark_verified_parameter_space_size, ...
    g_benchmark_verified_parameter_space_percentage] = print_parameter_box_size();
% total computation time
g_benchmark_total_time = toc;

switch (G_ANALYSIS_STRATEGY)
    case EOptions.ANALYSIS_ROVERGENE
        % just use the number of RoVerGeNe nodes
        g_benchmark_nodes_explored = g_benchmark_rovergene_exists_count;
    
    case EOptions.ANALYSIS_SPACEEX_ONLY
        % just use the number of SpaceEx nodes
        g_benchmark_nodes_explored = g_benchmark_hydentify_exists_count;
    
    case EOptions.ANALYSIS_EAGER
        % we have already incremented the count for each additional SpaceEx node
        g_benchmark_nodes_explored = ...
            g_benchmark_nodes_explored + g_benchmark_rovergene_exists_count;
    
    case EOptions.ANALYSIS_LAZY
        if (G_STORE_SEARCH_TREE)
            % the number of nodes can be read from the tree
            g_benchmark_nodes_explored = tree_size;
            % subtract nodes which were not explored
            for idx = tree_size : -1 : 1
                switch (tree{idx, EOptions.TREE_IDX_SKIP})
                    case {EOptions.SKIP_NORMAL, EOptions.SKIP_NORMAL_INFERRED}
                        % nothing to do;
                    
                    case {EOptions.SKIP_REDUNDANT, ...
                          EOptions.SKIP_REDUNDANT_INFERRED, ...
                          EOptions.SKIP_INCONSISTENT, ...
                          EOptions.SKIP_INCONSISTENT_INFERRED}
                        g_benchmark_nodes_explored = ...
                            g_benchmark_nodes_explored - 1;
                    
                    otherwise
                        assert(false, 'Unknown mode.');
                end
            end
        else
            % TODO implement this mode
        end
    
    otherwise
        assert(false, 'Unknown strategy.');
end

if (G_DEBUG_OUTPUT)
    fprintf('---------------\nBenchmark data:\n');
    fprintf('RoVerGeNe analyzed %d exists-KS in %s s.\n', ...
        g_benchmark_rovergene_exists_count, num2str(g_benchmark_rovergene_exists_time));
    fprintf('RoVerGeNe analyzed %d forall-KS in %s s.\n', ...
        g_benchmark_rovergene_forall_count, num2str(g_benchmark_rovergene_forall_time));
    fprintf('Hydentify analyzed %d exists-HA in %s s.\n', ...
        g_benchmark_hydentify_exists_count, num2str(g_benchmark_hydentify_exists_time));
    fprintf('Hydentify analyzed %d forall-HA in %s s.\n\n', ...
        g_benchmark_hydentify_forall_count, num2str(g_benchmark_hydentify_forall_time));
    fprintf('total analysis RoVerGeNe: %d KS in %s s.\n', ...
        g_benchmark_rovergene_exists_count + g_benchmark_rovergene_forall_count, ...
        num2str(g_benchmark_rovergene_exists_time + g_benchmark_rovergene_forall_time));
    fprintf('total analysis Hydentify: %d HA in %s s.\n', ...
        g_benchmark_hydentify_exists_count + g_benchmark_hydentify_forall_count, ...
        num2str(g_benchmark_hydentify_exists_time + g_benchmark_hydentify_forall_time));
    fprintf('parameter set size: %f = %f %% of %f\n', ...
        g_benchmark_verified_parameter_space_size, ...
        g_benchmark_verified_parameter_space_percentage, ...
        g_benchmark_total_parameter_space_size);
    if (g_benchmark_reachable_locations_spaceex_reachability_total > 0)
        fprintf('HyPy saved %d locations in total\n', ...
            g_benchmark_reachable_locations_spaceex_reachability_total);
    end
    fprintf('---------------\n');
end

%
% END OF FUNCTION
%

% ----- nested functions -----

% --- functions implementing different search strategies ---

function [ m_is_valid, is_descend_in_search_tree ] = strategy_rovergene(is_descend_in_search_tree)
% pure RoVerGeNe analysis

    % -- RoVerGeNe exists analysis --
    m_is_valid = rovergene_exists(transition_relation_existsP, ...
                                  transient_set_forallP);
    
    if (m_is_valid)
        % interpretation
        % parameter set is valid, stop here
        is_descend_in_search_tree = false;
    end
    
    % continue with forall analysis only if reasonable:
    % The only reason for forall analysis is to prune the search tree.
    if (is_descend_in_search_tree)
        % -- RoVerGeNe forall analysis --
        is_descend_in_search_tree = ...
            rovergene_forall(transition_relation_forallP, ...
                             transient_set_existsP);
    end
end

function [ m_is_valid, is_descend_in_search_tree, reachable_locations ] = strategy_spaceex_only(is_descend_in_search_tree, reachable_locations, spaceex_mode)
% pure SpaceEx analysis

    string = get_bit_string(b);
    
    % -- Hydentify exists analysis --
    [m_is_valid, reachable_locations] = hydentify_exists(Pb, string, ...
        transition_relation_existsP, reachable_locations, spaceex_mode);
    
    if (m_is_valid)
        % interpretation
        % parameter set is valid, stop here
        is_descend_in_search_tree = false;
    end
    
    % continue with forall analysis only if reasonable:
    % The only reason for forall analysis is to prune the search tree.
    if (is_descend_in_search_tree)
        % -- Hydentify forall analysis --
        is_descend_in_search_tree = hydentify_forall(Pb, string, ...
            transition_relation_forallP);
    end
end

function [ m_is_valid, is_descend_in_search_tree, reachable_locations ] = strategy_eager(is_descend_in_search_tree, reachable_locations, spaceex_mode)
% TODO describe

    % true iff parameter set is valid
    m_is_valid = false;

    % -- RoVerGeNe exists analysis --
    if ((CURRENT_MODE == EOptions.MODE_NORMAL) && ...
            rovergene_exists(transition_relation_existsP, transient_set_forallP))
        % interpretation
        % parameter set is valid, stop here
        m_is_valid = true;

    % -- Hydentify exists analysis --
    else
        string = get_bit_string(b);
        if (CURRENT_MODE == EOptions.MODE_SPACEEX_BOTH)
            g_benchmark_nodes_explored = g_benchmark_nodes_explored + 1;
        end
        [answer, reachable_locations] = hydentify_exists(Pb, string, ...
            transition_relation_existsP, reachable_locations, spaceex_mode);
        if (answer)
            % parameter set is valid, stop here
            m_is_valid = true;

            % output if a change towards RoverGene was found
            if (G_DEBUG_OUTPUT)
                fprintf('  <strong>Improvement by Hydentify in exists automaton!</strong>\n');
            end
        end
        clear answer;
    end

    % stop here when a valid result was found
    if (m_is_valid)
        is_descend_in_search_tree = false;
    end

    % continue with forall analysis only if reasonable:
    % The only reason for forall analysis is to prune the search tree.
    if (is_descend_in_search_tree)
        % -- RoVerGeNe forall analysis --
        if ((CURRENT_MODE == EOptions.MODE_NORMAL) && ...
                (rovergene_forall(transition_relation_forallP, ...
                                    transient_set_existsP)))
            % No intersection in RoverGene implies no intersection in
            % Hydentify, so there is no need to analyze in the other case.

            % -- Hydentify forall analysis --
        else
            answer = hydentify_forall(Pb, string, transition_relation_forallP);
            if (~ answer)
                % intersection found, no need to continue in this branch
                is_descend_in_search_tree = false;
            else
                % no intersection found, descend in the branch if reasonable,
                % but only with SpaceEx analysis
                if (G_DEBUG_OUTPUT && (CURRENT_MODE == EOptions.MODE_NORMAL))
                    % output if a change towards RoverGene was found
                    fprintf('  <strong>Improvement by Hydentify in forall automaton!</strong>\n');
                end
                
                CURRENT_MODE = EOptions.MODE_SPACEEX_BOTH;
            end
        end
    end
end

function [ is_valid, is_descend_in_search_tree, reachable_locations ] = strategy_lazy(is_descend_in_search_tree, reachable_locations, spaceex_mode)
% Do SpaceEx exists analysis only when
% 1) the RoVerGeNe exists analysis was successful (then
%    analyze the parent node)
% or
% 2) the RoVerGeNe exists analysis was not successful
%    and the current node is a leaf (then analyze the
%    current node).
% TODO describe

    % true iff parameter set is valid
    is_valid = false;
    
    switch (CURRENT_MODE)
        case EOptions.MODE_NORMAL
            % -- RoVerGeNe exists analysis --
            is_valid = rovergene_exists(transition_relation_existsP, ...
                                 transient_set_forallP);
            
            % push current or predecessor node to SpaceEx exists analysis stack
            if (is_valid || (B_LENGTH == CONSTRAINTS_LENGTH))
                push_to_spaceex_stack(b, is_valid, EOptions.MODE_SPACEEX_EX, ...
                                      true, node_idx);
                is_descend_in_search_tree = false;
                
            % -- RoVerGeNe forall analysis --
            elseif (~ rovergene_forall(transition_relation_forallP, ...
                                        transient_set_existsP))
                % RoVerGeNe cannot succeed in this branch, push current node to
                % SpaceEx analysis stack (both exists and forall)
                push_to_spaceex_stack(b, false, ...
                              EOptions.MODE_SPACEEX_BOTH_FIRST, true, node_idx);
                is_descend_in_search_tree = false;
            else
                assert(remember_node == -1, ...
                       'The remember node should be emptied.');
                remember_node = b;
                if (G_STORE_SEARCH_TREE)
                    remember_node_idx = node_idx;
                end
            end
        
        case EOptions.MODE_SPACEEX_EX
            % -- Hydentify exists analysis --
            [answer, reachable_locations] = hydentify_exists(Pb, ...
                get_bit_string(b), transition_relation_existsP, ...
                reachable_locations, spaceex_mode);
            if (answer)
            	is_valid = true;
                % push the predecessor to the SpaceEx exists analysis stack
                push_to_spaceex_stack(b, true, CURRENT_MODE, true, node_idx);
                if (G_DEBUG_OUTPUT)
                    % output if a change towards RoverGene was found
                    fprintf('  Improvement by Hydentify in exists automaton!\n');
                end
            end
            is_descend_in_search_tree = false;
        
        case {EOptions.MODE_SPACEEX_BOTH, EOptions.MODE_SPACEEX_BOTH_FIRST}
            string = get_bit_string(b);
            % -- Hydentify exists analysis --
            [answer, reachable_locations] = hydentify_exists(Pb, string, ...
                transition_relation_existsP, reachable_locations, spaceex_mode);
            if (answer)
                if (CURRENT_MODE == EOptions.MODE_SPACEEX_BOTH_FIRST)
                    % if this was the last node RoVerGeNe could explore and
                    % SpaceEx found a valid set, explore the predecessor
                    push_to_spaceex_stack(b, true, EOptions.MODE_SPACEEX_EX, ...
                                          true, node_idx);
                end
                
            	is_valid = true;
                is_descend_in_search_tree = false;
                if (G_DEBUG_OUTPUT)
                    % output if a change towards RoverGene was found
                    fprintf('  Improvement by Hydentify in exists automaton!\n');
                end
            % -- Hydentify forall analysis --
            elseif (is_descend_in_search_tree)
                answer = hydentify_forall(Pb, string, ...
                    transition_relation_forallP);
                if (~ answer)
                    is_descend_in_search_tree = false;
                else
                    % Hydentify tells us to continue search (but we can ignore
                    % RoVerGeNe from now on)
                    if (G_DEBUG_OUTPUT)
                        % output if a change towards RoverGene was found
                        fprintf('  Improvement by Hydentify in forall automaton!\n');
                    end
                end
            end
        
        otherwise
            assert(false, 'Unknown mode.');
    end
end

% --- data structure functions ---

function check_and_push_remember_node()
% in lazy strategy and RoVerGeNe mode, when a leaf is at the lowest possible
% level, push the transitive parent which was not redundant/inconsistent

    switch (G_ANALYSIS_STRATEGY)
        case {EOptions.ANALYSIS_ROVERGENE, EOptions.ANALYSIS_SPACEEX_ONLY, ...
              EOptions.ANALYSIS_EAGER}
            % nothing to do
        
        case EOptions.ANALYSIS_LAZY
            % push child nodes which are inconsistent
            if (B_LENGTH == CONSTRAINTS_LENGTH)
                if (G_DEBUG_OUTPUT)
                    fprintf('pushing transitive parent to SpaceEx stack:');
                    disp(remember_node);
                end
                
                assert(((~ G_STORE_SEARCH_TREE) || ...
                        (remember_node_idx ~= -1)), ...
                    'The tree node should have been stored as well.');
                push_to_spaceex_stack(remember_node, false, ...
                    EOptions.MODE_SPACEEX_EX, true, remember_node_idx);
                remember_node = -1;
                if (G_STORE_SEARCH_TREE)
                    remember_node_idx = -1;
                end
            end
        
        otherwise
            assert(false, 'Unknown strategy.');
    end
end

function push_to_normal_stack(b)
% push operation for normal stack

    stack_size = stack_size + 1;
    stack{stack_size} = b;
end

function [ b ] = pop_from_normal_stack()
% pop operation for normal stack

    b = stack{stack_size};
    stack_size = stack_size - 1;
end

function push_to_spaceex_stack(b, is_predecessor, mode, is_update_tree, node_idx)
% push operation for SpaceEx stack
% There is an option to add the predecessor node.

    if (is_predecessor)
        if (isempty(b))
            % no predecessor for root node
            return;
        else
            b = b(1 : length(b) - 1);
        end
        
        if (G_STORE_SEARCH_TREE)
            % Skip when the exists analysis in this node has already been
            % executed. This can happen when both child nodes are safe in
            % RoVerGeBe, but the parent was not.
            m_node_idx = tree{node_idx, EOptions.TREE_IDX_PARENT};
            if ((tree{m_node_idx, EOptions.TREE_IDX_DATA}...
                        (EOptions.DATA_IDX_SRE) ~= 0) || ...
                    (tree{m_node_idx, EOptions.TREE_IDX_STACK} < 0))
                if (G_DEBUG_OUTPUT)
                    fprintf('skipping node push %d to SpaceEx stack\n', ...
                        m_node_idx);
                end

                return;
            end
        end
    end
    
    stack_size_spaceex = stack_size_spaceex + 1;
    stack_spaceex{stack_size_spaceex, 1} = b;
    stack_spaceex{stack_size_spaceex, 2} = mode;
    
    if (G_STORE_SEARCH_TREE)
        if (is_update_tree)
            if (is_predecessor)
                set_node_stack_entry(tree{node_idx, EOptions.TREE_IDX_PARENT});
            else
                set_node_stack_entry(node_idx);
            end
        end
    end
    
    if (G_DEBUG_OUTPUT)
        if (is_predecessor)
            fprintf('pushing predecessor node');
        else
            fprintf('pushing current node');
        end
        if (G_STORE_SEARCH_TREE)
            if (is_predecessor)
                fprintf(' %d', tree{node_idx, EOptions.TREE_IDX_PARENT});
            else
                fprintf(' %d', node_idx);
            end
        end
        fprintf(' to SpaceEx stack, mode = %s\n', char(mode));
    end
end

function [ b, mode ] = pop_from_spaceex_stack()
% pop operation for SpaceEx stack

    switch (G_NODE_HEURISTICS)
        case EOptions.LAZY_SPACEEX_UNORDERED
            b = stack_spaceex{stack_size_spaceex, 1};
            mode = stack_spaceex{stack_size_spaceex, 2};
            stack_size_spaceex = stack_size_spaceex - 1;
        
        case {EOptions.LAZY_SPACEEX_CHILDREN, EOptions.LAZY_SPACEEX_PARENTS}
            m_idx = stack_size_spaceex;
            extreme = size(stack_spaceex{m_idx, 1}, 2);
            for m_i = 1 : (stack_size_spaceex - 1)
                current_size = size(stack_spaceex{m_i, 1}, 2);
                switch (G_NODE_HEURISTICS)
                    case EOptions.LAZY_SPACEEX_CHILDREN
                        swap = (current_size > extreme);
                    
                    case EOptions.LAZY_SPACEEX_PARENTS
                        swap = (current_size < extreme);
                    
                    otherwise
                        assert(false, 'unreachable code');
                end
                if (swap)
                    extreme = current_size;
                    m_idx = m_i;
                end
            end
            
            b = stack_spaceex{m_idx, 1};
            mode = stack_spaceex{m_idx, 2};
            
            if (m_idx ~= stack_size_spaceex)
                % swap with topmost element
                stack_spaceex{m_idx, 1} = stack_spaceex{stack_size_spaceex, 1};
                stack_spaceex{m_idx, 2} = stack_spaceex{stack_size_spaceex, 2};
                
                if (G_STORE_SEARCH_TREE)
                    % find node
                    m_node_idx = 1;
                    b_other = stack_spaceex{m_idx, 1};
                    for b_idx = b_other
                        if (b_idx)
                            m_node_idx = ...
                                tree{m_node_idx, EOptions.TREE_IDX_RIGHT};
                        else
                            m_node_idx = ...
                                tree{m_node_idx, EOptions.TREE_IDX_LEFT};
                        end
                    end
                    tree{m_node_idx, EOptions.TREE_IDX_STACK} = -m_idx;
                end
            end
            
            stack_size_spaceex = stack_size_spaceex - 1;
        
        otherwise
            assert(false, 'Unknown mode.');
    end
end

function remove_from_spaceex_stack(idx)
% removes an element at the given index from the stack

    % swap with topmost element
    stack_index = -tree{idx, EOptions.TREE_IDX_STACK};
    
    if (stack_index ~= stack_size_spaceex)
        % swap with topmost element
        stack_spaceex{stack_index, 1} = stack_spaceex{stack_size_spaceex, 1};
        stack_spaceex{stack_index, 2} = stack_spaceex{stack_size_spaceex, 2};

        % find node
        m_node_idx = 1;
        for b_idx = stack_spaceex{stack_index, 1}
            if (b_idx)
                m_node_idx = tree{m_node_idx, EOptions.TREE_IDX_RIGHT};
            else
                m_node_idx = tree{m_node_idx, EOptions.TREE_IDX_LEFT};
            end
        end

        % overwrite stack index
        tree{m_node_idx, EOptions.TREE_IDX_STACK} = -stack_index;
    end
    
    stack_size_spaceex = stack_size_spaceex - 1;
    tree{idx, EOptions.TREE_IDX_STACK} = 0;
    tree{idx, EOptions.TREE_IDX_DATA}(EOptions.DATA_IDX_SRE) = ...
        EOptions.DATA_NEG_INFERRED;
    
    if (G_DEBUG_OUTPUT)
        fprintf('removing node %d from analysis stack\n', idx);
    end
end

function remove_children_from_spaceex_stack(idx)
% remove child nodes from stack when SpaceEx exists analysis succeeded

    stack_index = tree{idx, EOptions.TREE_IDX_STACK};
    assert(stack_index <= 0, 'RoVerGeNe analysis should already be finished');
    if (stack_index < 0)
        remove_from_spaceex_stack(idx);
    end
    child = tree{idx, EOptions.TREE_IDX_LEFT};
    if (~ isempty(child))
        remove_children_from_spaceex_stack(child);
    end
    child = tree{idx, EOptions.TREE_IDX_RIGHT};
    if (~ isempty(child))
        remove_children_from_spaceex_stack(child);
    end
end

function remove_parents_from_spaceex_stack(idx)
% remove parent nodes from stack when SpaceEx exists analysis failed

    parent = tree{idx, EOptions.TREE_IDX_PARENT};
    while (~ isempty(parent))
        switch (tree{parent, EOptions.TREE_IDX_DATA}(EOptions.DATA_IDX_SRE))
            case {EOptions.DATA_NEG, EOptions.DATA_NEG_INFERRED}
                % when the entry is already negative, then we have already
                % propagated from here onwards
                assert(tree{parent, EOptions.TREE_IDX_STACK} == 0, ...
                    'The node should not be on the stack.');
                break;
            
            case EOptions.DATA_POS
                assert(false, ...
                    'The child node was unsafe, so the parent cannot be safe.');
            
            case 0
                % skip over this node
            
            otherwise
                assert(false, 'unknown result');
        end
        
        stack_index = tree{parent, EOptions.TREE_IDX_STACK};
        assert(stack_index <= 0, ...
            'RoVerGeNe analysis should already be finished');
        if (stack_index < 0)
            remove_from_spaceex_stack(parent);
        end
        % set inferred result to avoid future computations
        tree{parent, EOptions.TREE_IDX_DATA}(EOptions.DATA_IDX_SRE) = ...
                EOptions.DATA_NEG_INFERRED;
        
        parent = tree{parent, EOptions.TREE_IDX_PARENT};
    end
end

function add_tree_node(stack_index, parent, b, reachable_locations)
% adds a node to the search tree

    tree_size = tree_size + 1;
    tree{tree_size, EOptions.TREE_IDX_STACK} = stack_index;
    tree{tree_size, EOptions.TREE_IDX_SKIP} = 0;
    tree{tree_size, EOptions.TREE_IDX_DATA} = zeros(4, 1);
    tree{tree_size, EOptions.TREE_IDX_REACHABLE_LOCATIONS} = reachable_locations;
    tree{tree_size, EOptions.TREE_IDX_BOOLVEC} = b; % TODO only for debugging
    
    % update parent and child pointers (negative iff rhs child)
    if (parent < 0)
        parent = -parent;
        tree{parent, EOptions.TREE_IDX_RIGHT} = tree_size;
        tree{tree_size, EOptions.TREE_IDX_PARENT} = parent;
    elseif (parent > 0)
        tree{parent, EOptions.TREE_IDX_LEFT} = tree_size;
        tree{tree_size, EOptions.TREE_IDX_PARENT} = parent;
    end
end

function set_node_stack_entry(m_idx)
% updates the stack entry of a search tree node to current topmost SpaceEx node

    assert((tree{m_idx, EOptions.TREE_IDX_STACK} == 0), ...
        'node appears in both stacks');
    tree{m_idx, EOptions.TREE_IDX_STACK} = -stack_size_spaceex;
end

% --- helper functions to call RoVerGeNe/Hydentify exists/forall analysis ---

function [ answer ] = rovergene_exists(transition_relation_existsP, transient_set_forallP)
    if (G_DEBUG_OUTPUT)
        fprintf('KS exists check started: |%g s|', toc);
    end
    timer = tic;
    answer = export_and_check_in_rovergene(1, transition_relation_existsP, ...
                                transient_set_forallP);
    g_benchmark_rovergene_exists_time = ...
        g_benchmark_rovergene_exists_time + toc(timer);
    g_benchmark_rovergene_exists_count = g_benchmark_rovergene_exists_count + 1;
    if (G_DEBUG_OUTPUT)
        fprintf(' -- finished with answer %d: |%g s|\n', answer, toc);
    end
    if (G_STORE_SEARCH_TREE)
        if (answer)
            tree{node_idx, EOptions.TREE_IDX_DATA}(EOptions.DATA_IDX_RGE) = ...
                EOptions.DATA_POS;
        else
            tree{node_idx, EOptions.TREE_IDX_DATA}(EOptions.DATA_IDX_RGE) = ...
                EOptions.DATA_NEG;
        end
    end
end

function [ answer ] = rovergene_forall(transition_relation_forallP, transient_set_existsP)
    if (G_DEBUG_OUTPUT)
        fprintf('KS forall check started: |%g s|', toc);
    end
    timer = tic;
    answer = export_and_check_in_rovergene(0, transition_relation_forallP, ...
                                transient_set_existsP);
    g_benchmark_rovergene_forall_time = ...
        g_benchmark_rovergene_forall_time + toc(timer);
    g_benchmark_rovergene_forall_count = g_benchmark_rovergene_forall_count + 1;
    if (G_DEBUG_OUTPUT)
        fprintf(' -- finished with answer %d: |%g s|\n', answer, toc);
    end
    if (G_STORE_SEARCH_TREE)
        if (answer)
            tree{node_idx, EOptions.TREE_IDX_DATA}(EOptions.DATA_IDX_RGA) = ...
                EOptions.DATA_POS;
        else
            tree{node_idx, EOptions.TREE_IDX_DATA}(EOptions.DATA_IDX_RGA) = ...
                EOptions.DATA_NEG;
        end
    end
end

function [ answer, reachable_locations ] = hydentify_exists(Pb, string, transition_relation, reachable_locations, spaceex_mode)
    if (G_DEBUG_OUTPUT)
        fprintf('LHA exists check started: |%g s|', toc);
    end
    timer = tic;
    [answer, reachable_locations] = export_and_check_in_spaceex(Pb, 1, string, ...
        transition_relation, reachable_locations, spaceex_mode);
    g_benchmark_hydentify_exists_time = ...
        g_benchmark_hydentify_exists_time + toc(timer);
    g_benchmark_hydentify_exists_count = ...
        g_benchmark_hydentify_exists_count + 1;
    if (G_DEBUG_OUTPUT)
        fprintf(' -- finished with answer %d: |%g s|\n', answer, toc);
    end

    % interpretation: "unknown" result due to SpaceEx
    if (answer == 2)
        if (G_DEBUG_OUTPUT)
            fprintf('\n!!!\n!!!\n  problem with exists automaton\n');
        end
        answer = 0;
    end
    if (G_STORE_SEARCH_TREE)
        if (answer)
            tree{node_idx, EOptions.TREE_IDX_DATA}(EOptions.DATA_IDX_SRE) = ...
                EOptions.DATA_POS;
            
            % remove child nodes from stack
            remove_children_from_spaceex_stack(node_idx);
        else
            tree{node_idx, EOptions.TREE_IDX_DATA}(EOptions.DATA_IDX_SRE) = ...
                EOptions.DATA_NEG;
            
            % remove parent nodes from stack
            remove_parents_from_spaceex_stack(node_idx);
        end
    end
end

function [ answer ] = hydentify_forall(Pb, string, transition_relation)
    % only analyze forall automaton when the RoverGene forall-TS found
    % an intersection
    if (G_DEBUG_OUTPUT)
        fprintf('LHA forall check started: |%g s|', toc);
    end
    timer = tic;
    [answer, ~] = export_and_check_in_spaceex(Pb, 0, string, ...
        transition_relation, [], EOptions.SPACEEX_EARLY_TERMINATION);
    g_benchmark_hydentify_forall_time = ...
        g_benchmark_hydentify_forall_time + toc(timer);
    g_benchmark_hydentify_forall_count = ...
        g_benchmark_hydentify_forall_count + 1;
    if (G_DEBUG_OUTPUT)
        fprintf(' -- finished with answer %d: |%g s|\n', answer, toc);
    end

    % interpretation: "unknown" result due to SpaceEx
    if (answer == 2)
        if (G_DEBUG_OUTPUT)
            fprintf('\n!!!\n!!!\n  problem with forall automaton\n');
        end
        answer = 1;
    end
    if (G_STORE_SEARCH_TREE)
        if (answer)
            tree{node_idx, EOptions.TREE_IDX_DATA}(EOptions.DATA_IDX_SRA) = ...
                EOptions.DATA_POS;
        else
            tree{node_idx, EOptions.TREE_IDX_DATA}(EOptions.DATA_IDX_SRA) = ...
                EOptions.DATA_NEG;
        end
    end
end

end % end of main function

% ----- other functions -----

function [ string ] = get_bit_string(b)
% convert bit array into string
    string = '[';
    for idx = b
        string = strcat(string, num2str(idx));
    end
    string = strcat(string, ']');
end

function [ total_size, found_size, percentage ] = print_parameter_box_size()
% prints the resulting parameter set size for simple boxes

    global unknown_parameter;
    global unknown_parameter_nb;
    global production_rate_parameter;
    global degradation_rate_parameter;
    
    found_size = 0;
    
    if (unknown_parameter_nb == 0)
        total_size = 0;
        percentage = 0;
        return;
    end
    
    % complete valid parameter set size
    global result;
    certain = result.valid_parameter_set.certain;
    for i = 1 : length(certain)
        polytope = certain(i);
        [H, K] = double(polytope);
        box = zeros(size(H, 2), 2);
        for j = 1 : length(K)
            row = H(j, :);
            for k = 1 : length(row)
                elem = row(k);
                if (elem ~= 0)
                    if (elem < 0)
                        box(k, 1) = -K(j);
                    else
                        box(k, 2) = K(j);
                    end
                    break;
                end
            end
        end
        
        box_size = 1;
        for j = 1 : size(box, 1)
            box_size = box_size * (box(j, 2) - box(j, 1));
        end
        found_size = found_size + box_size;
    end
    
    % compute total parameter set size
    total_size = 1;
    for i = 1 : unknown_parameter_nb
        indices = unknown_parameter{i};
        if (indices(1) == 1)
            parameter = production_rate_parameter{indices(2)}{indices(3)};
        else
            assert((indices(1) == -1), 'Unknown encoding of parameters.');
            parameter = degradation_rate_parameter{indices(2)}{indices(3)};
        end
        assert((size(parameter, 2) == 2), 'A parameter should be an interval.');
        total_size = total_size * (parameter(2) - parameter(1));
    end
    
    % compute percentage
    percentage = (found_size / total_size * 100);
end