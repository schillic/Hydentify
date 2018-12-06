function [ answer, reachable_locations ] = export_and_check_in_spaceex(Pb, EXISTS, FILENAME_SUFFIX, TRANSITION_RELATION_SPARSE, reachable_locations, SPACEEX_MODE)
%export_and_check_in_spaceex - export LHA and check in SpaceEx
% -------------------------------------------------------------------------
% Authors: Ezio Bartocci
%          SUNY at Stony Brook
%          Stony Brook, NY, USA
% email:   eziobart@ams.sunysb.edu
% Website: http://www.eziobartocci.com
%
%          Gregory Batt
%          Boston University, 
%          Brookline, MA, USA
% email:   batt@bu.edu
% Website: http://iasi.bu.edu/~batt/
%
%          Radu Grosu
%          SUNY at Stony Brook 
%          Stony Brook, NY, USA
% email:   grosu@cs.sunysb.edu
% Website: http://www.cs.sunysb.edu/~grosu/
%
% Author:  Christian Schilling
%          University of Freiburg, Germany
% email:   schillic at informatik.uni-freiburg.de
% Website: http://swt.informatik.uni-freiburg.de/staff/christian_schilling
%
% September 2015
% -------------------------------------------------------------------------

    global spaceex_file_name;
    global G_USE_OLD_RESULTS;
    global G_SKIP_EMPTINESS_TEST_FOR_FORALL_LHA;
    global G_SAMPLING_BOUND;
    global G_SAMPLING_MAX;
    global G_SAMPLING_SHRUNK;
    global G_SAMPLING_EMPTY;
    global G_SAMPLING_MAX_TOTAL;
    global G_SAMPLING_SHRUNK_TOTAL;
    global G_SAMPLING_EMPTY_TOTAL;
    G_SAMPLING_MAX = 0;
    G_SAMPLING_SHRUNK = 0;
    G_SAMPLING_EMPTY = 0;

    % options to skip export and/or check (e.g., reuse old results)
    export = true;
    check = true;

    % test emptiness for forall automaton (would lead to an error)
    if ((~ G_SKIP_EMPTINESS_TEST_FOR_FORALL_LHA) && (~ EXISTS))
        if (~ wrap_isempty(Pb))
            % For an empty parameter space we say that the forall-automaton
            % has no intersection.
            answer = 0;

            output('\nparameter space is empty, skipping forall automaton!\n\n');
            return;
        end
    end

    % file names
    if (EXISTS)
        FILENAME_RAW = strcat(spaceex_file_name.existsP, FILENAME_SUFFIX);
    else
        FILENAME_RAW = strcat(spaceex_file_name.forallP, FILENAME_SUFFIX);
    end
    PREFIX = './exports/';
    FILEPATHNAME_LHA = strcat(PREFIX, FILENAME_RAW, '.xml');
    FILEPATHNAME_RESULT = strcat(PREFIX, FILENAME_RAW, '-result.txt');

    % search for former output on disk
    if (G_USE_OLD_RESULTS)
        if (exist(FILEPATHNAME_RESULT, 'file') == 2)
            output('skipping check, using old result\n');
            export = false;
            check = false;
            fid = fopen(FILEPATHNAME_RESULT, 'r');
            answer = str2double(fgetl(fid));
            fclose(fid);
        elseif (exist(FILEPATHNAME_LHA, 'file') == 2)
            export = false;
        end
    end

    % export hybrid automaton file in SpaceEx format
    if (export)
        [initial_names, forbidden_names] = ...
            spaceex_export_LHA(Pb, EXISTS, TRANSITION_RELATION_SPARSE, ...
                               FILEPATHNAME_LHA, reachable_locations);
        if ((~ EXISTS) && (G_SAMPLING_BOUND > 0))
            fprintf('  <strong>needed %d-convergence for the whole automaton\n  shrunk %d location flows, %d of which are now empty</strong>\n', ...
                G_SAMPLING_MAX, G_SAMPLING_SHRUNK, G_SAMPLING_EMPTY);
            G_SAMPLING_MAX_TOTAL = max(G_SAMPLING_MAX, G_SAMPLING_MAX_TOTAL);
            G_SAMPLING_SHRUNK_TOTAL = max(G_SAMPLING_SHRUNK, G_SAMPLING_SHRUNK_TOTAL);
            G_SAMPLING_EMPTY_TOTAL = max(G_SAMPLING_EMPTY, G_SAMPLING_EMPTY_TOTAL);
        end
    else
        output('\nskipping export, using old automaton file\n');
    end
    
    % check whether any forbidden location is reachable in the discrete graph
    % NOTE: This should only occur when not run with RoVerGeNe, which also does
    % graph reachability for reachability properties.
    if (exist('forbidden_names', 'var'))
        if ((size(forbidden_names, 1) == 1) && isempty(forbidden_names{1}))
            output('skipping check, forbidden states are not reachable in the discrete structure\n');
            answer = 1;
            check = false;
        end
    end

    % check whether property is true using SpaceEx
    if (check)
        switch (SPACEEX_MODE)
            case EOptions.SPACEEX_EARLY_TERMINATION
                if (export)
                    answer = spaceex_check_property(FILENAME_RAW, ...
                        FILEPATHNAME_RESULT, 0, initial_names, forbidden_names);
                else
                    answer = spaceex_check_property(FILENAME_RAW, ...
                        FILEPATHNAME_RESULT, 1);
                end
                
            case EOptions.SPACEEX_HYPY
                if (export)
                    [answer, reachable_locations] = hypy_check_property(...
                        FILENAME_RAW, FILEPATHNAME_RESULT, initial_names, ...
                        reachable_locations);
                else
                    [answer, reachable_locations] = hypy_check_property(...
                        FILENAME_RAW, FILEPATHNAME_RESULT, [], ...
                        reachable_locations);
                end
                
            otherwise
                assert(false, 'Unknown mode.');
        end
    end
end


function [ initial_names, forbidden_names, id_offset_map ] = spaceex_export_LHA(Pb, EXISTS, ...
        TRANSITION_RELATION_SPARSE, FILEPATHNAME_LHA, reachable_locations)
% exports SpaceEx model file to folder "exports"

	global box_nb_i;
	global polytope_lib;
    global variable_nb;
    global G_INITIAL_LOCS;
    global G_FORBIDDEN_LOCS;
    global G_DEBUG_OUTPUT;
    
    % stimulus related variables
    global stimulus_mode;
    global stimulus_wrapper;
    switch (stimulus_mode)
        case {1, 2}
            IS_IGNORE_STIM = 1;
            STIM_DIM = stimulus_wrapper{1};
        otherwise
            IS_IGNORE_STIM = 0;
            STIM_DIM = 0;
    end
    
    % numbers of locations/variables
    NUM_LOC = prod(box_nb_i, 2);
    
    % convert sparse matrix into a simple pair mapping matrix
    [TRANS_REL(:, 1), TRANS_REL(:, 2)] = find(TRANSITION_RELATION_SPARSE);
    
    % detect unreachable locations in discrete structure to skip them
    if (~ isempty(G_INITIAL_LOCS))
        REACHED = graph_reachability(TRANS_REL, G_INITIAL_LOCS, NUM_LOC, ...
            G_DEBUG_OUTPUT, reachable_locations);
        initial_names = cell(size(G_INITIAL_LOCS, 2), 1);
        initial_location_names_ctr = 1;
        id_offset_map = zeros(NUM_LOC, 1);
        id_offset = 0;
    else
        REACHED = [];
        initial_names = cell(0);
        initial_location_names_ctr = 0;
        id_offset_map = -1;
    end
    if (~ isempty(G_FORBIDDEN_LOCS))
        if (isempty(REACHED))
            forbidden_names = cell(size(G_FORBIDDEN_LOCS, 2), 1);
            forbidden_location_idx = 1;
        else
            forbidden_location_idx = 0;
            forbidden_reached = 0;
            for i = 1 : size(G_FORBIDDEN_LOCS, 2)
                if (REACHED(G_FORBIDDEN_LOCS(1, i)) == 1)
                    if (forbidden_location_idx == 0)
                        forbidden_location_idx = i;
                    end
                    forbidden_reached = forbidden_reached + 1;
                end
            end
            if (forbidden_reached == 0)
                forbidden_names = cell(1);
                forbidden_names{1} = [];
                return;
            end
            forbidden_names = cell(forbidden_reached, 1);
        end
        forbidden_location_names_ctr = 1;
    else
        forbidden_names = cell(0);
        forbidden_location_names_ctr = 0;
    end
    
    % epsilon strings for fast (de-)activation
    ADD_EPSILONS = false;
    if (ADD_EPSILONS)
        EPS_PLUS = ' + 0.0001'; % epsilon for state variables
        EPS_MINUS = ' - 0.0001'; % epsilon for state variables
        EPS_TIME = ' + 0.0001'; % epsilon for time variable
    else
        EPS_PLUS = '';
        EPS_MINUS = '';
        EPS_TIME = '';
    end
	
    % output model file; open document and component
	FILE = fopen(FILEPATHNAME_LHA, 'w');
	spaceex_tag_open(FILE);
	component_tag_open(FILE, 'pma');
	
	% write parameters and create flow string
    % (which is constant, the actual flow is encoded in the invariant)
	FLOW_STRING = write_spaceex_params_and_flow(FILE, TRANS_REL, ...
        variable_nb, IS_IGNORE_STIM, STIM_DIM);
    
    % indices array to reason about current location
    indices = ones(variable_nb, 1);
    increment = 0;
    
    % write locations
    for i_box = 1 : NUM_LOC
        if (~ isempty(REACHED))
            % skip locations which are unreachable in the discrete system
            % For this we need to adapt the location IDs.
            if (REACHED(i_box) ~= 1)
                increment = increment + 1;
                id_offset = id_offset + 1;
                continue;
            else
                BOX_ID = (i_box - id_offset);
                id_offset_map(i_box) = BOX_ID;
            end
        else
            BOX_ID = i_box;
        end
        
        % compute current indices array
        [indices, increment] = get_next_index(indices, increment, variable_nb);
        
        % write location header
        box_name = write_location_header(FILE, BOX_ID, indices, variable_nb);
        
        % add to initial/forbidden locations
        if (initial_location_names_ctr > 0)
            if (G_INITIAL_LOCS(1, initial_location_names_ctr) == i_box)
                initial_names{initial_location_names_ctr} = box_name;
                if (initial_location_names_ctr == ...
                        size(initial_names, 1))
                    initial_location_names_ctr = 0;
                else
                    initial_location_names_ctr = initial_location_names_ctr + 1;
                end
            end
        end
        if (forbidden_location_names_ctr > 0)
            if (G_FORBIDDEN_LOCS(1, forbidden_location_idx) == i_box)
                forbidden_names{forbidden_location_names_ctr} = box_name;
                if (forbidden_location_names_ctr == size(forbidden_names, 1))
                    forbidden_location_names_ctr = 0;
                else
                    forbidden_location_names_ctr = ...
                        forbidden_location_names_ctr + 1;
                    
                    if (~ isempty(REACHED))
                        % find the next non-empty forbidden state
                        while (true)
                            forbidden_location_idx = ...
                                forbidden_location_idx + 1;

                            if (forbidden_location_idx > ...
                                    size(G_FORBIDDEN_LOCS, 2))
                                forbidden_location_names_ctr = 0;
                                break;
                            end

                            if (REACHED(G_FORBIDDEN_LOCS(1, forbidden_location_idx)) == 1)
                                break;
                            end
                        end
                    end
                end
            end
        end
        
        % compute polytope
        P = get_polytope(EXISTS, variable_nb, indices, Pb);
        
        % compute invariant
        inv_string = get_invariant_string(EXISTS, variable_nb, indices, P, ...
            IS_IGNORE_STIM, STIM_DIM, EPS_PLUS, EPS_MINUS, EPS_TIME);
        
        % write flow and invariant
        flow_tags(FILE, FLOW_STRING);
        invariant_tags(FILE, inv_string);
        
        % finalize location
        location_tag_close( FILE);

        % back to MPT
        polytope_lib = 'mpt';
    end
    
    % write transitions
    write_transitions(FILE, TRANS_REL, id_offset_map);
    
    % finalize component and system
    component_tag_close(FILE);
	spaceex_tag_close(FILE);
	fclose(FILE);
end


function [ flow ] = write_spaceex_params_and_flow(FILE, TRANS_REL, NUM_VAR, IS_IGNORE_STIM, STIM_DIM)
% declares the variables etc. in the SpaceEx model and constructs the flow
% string

    flow = '';

    for i = 1 : NUM_VAR
        % ignore stimulus dimension
        if (IS_IGNORE_STIM && i == STIM_DIM)
            continue;
        end

        numstring = num2str(i);

        % variables xi
        param_tag(FILE, strcat('x', numstring), 'real', 'false', '1', '1', ...
                  'any');

        % variables fi
        param_tag_controlled(FILE, strcat('f', numstring), 'real', 'true', ...
                             '1', '1', 'any', 'false');

        % flows
        flow = strcat(flow, ' x', numstring, ''' == f', numstring, ' &amp;');
    end

    global G_ADD_TIME;
    if (G_ADD_TIME)
        % add time variable
        param_tag(FILE, 't', 'real', 'true', '1', '1', 'any');
        flow = strcat(flow, ' t''==1');
    else
        % remove last ampersand (substring " &amp;")
        flow = flow(1 : length(flow) - 6);
    end

    % add max_time
%     param_tag ( fid, 'max_time', 'real', 'false', '1', '1', 'const');

    % add transition epsilons if enabled
    global G_TRANSITION_EPSILONS;
    if (G_TRANSITION_EPSILONS > 0)
        param_tag ( FILE, 'eps', 'real', 'false', '1', '1', 'const');
        param_tag ( FILE, 't_eps', 'real', 'true', '1', '1', 'any');
        flow = strcat(flow, ' &amp; t_eps''==1');
    end

    % labels
    for i = 1 : size(TRANS_REL, 1)
        src = TRANS_REL(i, 1);
        dest = TRANS_REL(i, 2);
        if (src ~= dest) % skip the self-loops
            param_tag_label(FILE, sprintf('T%d_%d', src, dest), 'true');
        end
    end
end


function [ indices, increment ] = get_next_index(indices, increment, NUM_VAR)
% increments the indices array by the given number

    global box_nb_i;
    
    for inc = 1 : increment
        for i = 1 : NUM_VAR
            if (indices(i) < box_nb_i(i))
                indices(i) = indices(i) + 1;
                break;
            else
                indices(i) = 1;
            end
        end
    end
    increment = 1;
end


function [ box_name ] = write_location_header(FILE, box_id, indices, NUM_VAR)
% writes the location header and generates the indices array

    global box_nb_i;
    
    % construct location name and generate indices array
    box_name = 'Box';
    for idx = 1 : NUM_VAR
        box_name = strcat(box_name, '_', num2str(indices(idx)));
    end
    
    % TODO should we have a mode where all locations have constant positions?
    % TODO is there a simpler way to place the locations?
    
    % standard location sizes
	BOX_SIZE  = 256;
	BOX_SPACE = 64;
    SIZE_PLUS_SPACE = BOX_SIZE + BOX_SPACE;
    
    % X position = dimension 1
    X_POS = BOX_SPACE + (BOX_SIZE / 2) + SIZE_PLUS_SPACE * (indices(1) - 1);
    % Y position = all other dimensions
    y_pos = - (BOX_SIZE / 2);
    if (NUM_VAR > 1)
        product = 1;
        for dim = 2 : NUM_VAR
            y_pos = y_pos - (SIZE_PLUS_SPACE * (indices(dim) - 1) * product);
            product = product * box_nb_i(dim);
        end
        y_pos = y_pos + (SIZE_PLUS_SPACE * product);
    end
    location_tag_open(FILE, box_id, box_name, X_POS, y_pos, BOX_SIZE, BOX_SIZE);
end


function [ P ] = get_polytope(EXISTS, NUM_VAR, indices, Pb)
% computes the polytope for the given parameter space
% The computation and result differ for the exists and the forall automaton.

    global partition;
    global polytope_lib;
    
    % compute the vertices of the polytope of Pspace
    VERTICES_P = wrap_extreme(Pb);

    % use PPL for the rest
    polytope_lib = 'pplmex';

    % compute the number of vertices in Pspace
    NUM_VERTICES_P = size(VERTICES_P, 1);
    
    % number of vertices of a state space box = 2^{num of the dimensions of box}
    NUM_VERTICES_S = 2^NUM_VAR;
    
    % calculate the state space box
    vertices_s = zeros(NUM_VERTICES_S, NUM_VAR);
    for dim = 1 : NUM_VAR
        for idx_s = 1 : NUM_VERTICES_S
            vertices_s(idx_s, dim) = partition{dim}...
                (indices(dim) + mod(floor((idx_s - 1) / 2^(dim - 1)), 2));
        end
    end

    % compute P with previously computed data
    if (EXISTS)
        P = get_polytope_exists(NUM_VAR, NUM_VERTICES_S, NUM_VERTICES_P, ...
                                vertices_s, VERTICES_P);
    else
        P = get_polytope_forall(NUM_VAR, NUM_VERTICES_S, NUM_VERTICES_P, ...
                                vertices_s, VERTICES_P, Pb);
    end
end


function [ P ] = get_polytope_exists(NUM_VAR, NUM_VERTICES_S, NUM_VERTICES_P, VERTICES_S, VERTICES_P)
% computes the polytope for the given parameter space for an exists automaton

    if (NUM_VERTICES_P == 0)
        % empty parameter space
        START = 0;
        vertices_all = zeros(NUM_VERTICES_S, NUM_VAR);
    else
        % non-empty parameter space
        START = 1;
        vertices_all = zeros(NUM_VERTICES_S * NUM_VERTICES_P, NUM_VAR);
    end
    
    ind = 0;
    for idx_p = START : NUM_VERTICES_P
        param = VERTICES_P(idx_p, :);
        for idx_s = 1 : NUM_VERTICES_S
            f_dim = get_polytope_evaluation_loop(NUM_VAR, ...
                                VERTICES_S(idx_s, :), param);
            ind = ind + 1;
            vertices_all(ind, :) = f_dim;
        end
    end
	
	P = wrap_polytope(vertices_all);
end


function [ intersections ] = get_polytope_forall(NUM_VAR, NUM_VERTICES_S, NUM_VERTICES_P, VERTICES_S, VERTICES_P, Pb)
% computes the polytope for the given parameter space for a forall automaton

	vertices_all = zeros(NUM_VERTICES_S, NUM_VAR);
    is_empty = false;
	for idx_p = 1 : NUM_VERTICES_P
        param = VERTICES_P(idx_p, :);
        for idx_s = 1 : NUM_VERTICES_S
            f_dim = get_polytope_evaluation_loop(NUM_VAR, ...
                                VERTICES_S(idx_s, :), param);
            vertices_all(idx_s, :) = f_dim;
        end
		% intersect over all boxes
        P = wrap_polytope(vertices_all);
		if (idx_p == 1)
			intersections = P;
		else
			intersections = P & intersections;
            
            % stop if polytope is empty
            if (wrap_isempty(intersections))
                is_empty = true;
                break;
            end
		end
	end
    
    % add sampling for tighter approximation
    global G_SAMPLING_BOUND;
    if ((G_SAMPLING_BOUND > 0) && (~ is_empty))
        intersections = sample(intersections, Pb, NUM_VAR, NUM_VERTICES_S, ...
            VERTICES_S);
    end
end


function [ f_dim ] = get_polytope_evaluation_loop(NUM_VAR, vertex, param)
% common part for computing the polytope for the given parameter space

    f_dim = zeros(NUM_VAR, 1);
    % for each dimension there is a function to evaluate
    for var = 1 : NUM_VAR
        % Check whether the vertex has the lower or the upper
        % value in this box (needed for the step function mode).
        % NOTE: currently not used
%         if (init_box_v(1, var) == vertex(var))
%             % init_box_v(ind, 0) has the lower entry for each dimension
%             isLower = 1;
%         else
%             assert(init_box_v(size(init_box_v, 1), var) == vertex(var), ...
%                 'The vertex must be a corner of the box.');
%             isLower = 0;
%         end
        [f_dim(var), ~] = evaluate_f(vertex, var, 1, param);
    end
end


function [ intersections ] = sample(intersections, Pb, NUM_VAR, NUM_VERTICES_S, INIT_VERTICES)
% add intersection with sampled points to be more precise in forall flow

    global G_SAMPLING_BOUND;
    global G_SAMPLING_MAX;
    global G_SAMPLING_SHRUNK;
    global G_SAMPLING_EMPTY;
    noLongestTry = 0;
    noLongestBoundNeeded = 0;
    noSamples = 0;
    isShrunk = 0;
    vertices_all = zeros(NUM_VERTICES_S, NUM_VAR);
    while (noLongestTry < G_SAMPLING_BOUND)
        % use MPT's random function
        for idx_s = 1 : NUM_VERTICES_S
            f_dim = get_polytope_evaluation_loop(NUM_VAR, ...
                                        INIT_VERTICES(idx_s, :), randpoint(Pb));
            vertices_all(idx_s, :) = f_dim;
        end
        
        P = wrap_polytope(vertices_all);
        intersections_new = P & intersections;
        if (eq(intersections_new, intersections))
            noLongestTry = noLongestTry + 1;
        else
            isShrunk = 1;
            intersections = intersections_new;
            noLongestBoundNeeded = noLongestTry;
            noLongestTry = 0;
        end
        noSamples = noSamples + 1;
    end
    
    if (isShrunk)
        if (wrap_isempty(intersections))
            G_SAMPLING_EMPTY = G_SAMPLING_EMPTY + 1;
        end
        G_SAMPLING_SHRUNK = G_SAMPLING_SHRUNK + 1;
    end
    
    G_SAMPLING_MAX = max(G_SAMPLING_MAX, noLongestBoundNeeded);
%     fprintf('total samples: %d, needed %d-convergence\n', noSamples, ...
%         noLongestBoundNeeded);
end


function [ inv_string ] = get_invariant_string(exist, NUM_VAR, indices, ...
        P, ignore_stimulus, stimulus_dimension, epsilon_plus, epsilon_minus, ...
        epsilon_time)
% constructs the invariant string

    global partition;
    global box_nb_i;
    global time_constraint;
    global stimulus_mode;
    global stimulus_wrapper;
    global polytope_scale;
    global G_EXACT_FLOW; % parameter switch for exact computation
    global G_ADD_TIME;
    
    if (G_EXACT_FLOW)
        % exact constraints: polynomial(X') <= f <= polynomial(X'), where X' is
        % all variables without x.
        
        % The factors are only correct for the vector. That is why below we
        % multiply again.
        [H, K] = wrap_hk(P);
        H = H * polytope_scale;
        size_w = size(H, 2);
        
        if (exist || (size_w > 0))
            inv1 = '';
        else
            % forall automaton: intersection is empty -> invariant is 'false'
            inv_string = '1==0';
            return;
        end
        
        if size_w > 0
            is_printed_line = false;
            for line = 1 : size(K, 1)
                % ignore the constraints for the stimulus (row)
                if ((ignore_stimulus) && (H(line, stimulus_dimension) ~= 0))
                    continue;
                end
                
                lhs = '';
                is_printed_column = false;
                for column = 1 : size(H, 2)
                    % ignore stimulus dimension
                    if ((ignore_stimulus) && (column == stimulus_dimension))
                        % ignore the stimulus dimension (column)
                        continue;
                    end
                    
                    if (is_printed_column)
                        lhs = sprintf('%s + ', lhs);
                    else
                        is_printed_column = true;
                    end
                    lhs = sprintf('%s%d * f%d', lhs, H(line, column), column);
                end
                
                if (is_printed_line)
                    inv1 = sprintf('%s &amp;\n', inv1);
                else
                    is_printed_line = true;
                end
                inv1 = sprintf('%s %s &lt;= %d', inv1, lhs, K(line, 1));
            end
        end
    else
        % lower bound <= f <= upper bound
        
        % We manage a cube that contains our dynamics
        % In Phaver we can express the dynamics in terms of
        % fx_min <= fx <= fx_max
        dyn_box_points = wrap_extreme(P);
        size_w = size(dyn_box_points, 2);
        
        if (exist || (size_w > 0))
            inv1 = '';
        else
            % forall automaton: intersection is empty -> invariant is 'false'
            inv_string = '1==0';
            return;
        end
        
        past_first = 0;
        for d = 1 : size_w
            % ignore stimulus dimension
            if (ignore_stimulus && d == stimulus_dimension)
                continue;
            end
            
            if (past_first)
                inv1 = sprintf('%s &amp;\n', inv1);
            else
                past_first = 1;
            end
            
            inv1 = sprintf('%s %s &lt;= f%d &lt;= %s', ...
                inv1, num2str_wrap(min(dyn_box_points(:, d))), d, ...
                num2str_wrap(max(dyn_box_points(:, d))));
        end
    end
    
    % independent part
    inv2 = '';
    for m_i = 1 : NUM_VAR
        % ignore stimulus dimension
        if (ignore_stimulus && m_i == stimulus_dimension)
            % ignore stimulus dimension if enabled
            continue;
        end
        index = indices(m_i);
        partition_i = partition{m_i};
        if (m_i == stimulus_dimension)
            % no epsilons for stimulus dimension
            inv2 = sprintf('%s %s &lt;= x%d &lt;= %s &amp;\n', ...
                inv2, num2str_wrap(partition_i(index)), m_i, ...
                num2str_wrap(partition_i(index + 1)));
        elseif (index == 1)
            % no lower epsilons for leftmost boxes
            inv2 = sprintf('%s %s &lt;= x%d &lt;= %s %s &amp;\n', ...
                inv2, num2str_wrap(partition_i(index)), m_i, ...
                num2str_wrap(partition_i(index + 1)), epsilon_plus);
        elseif (index == box_nb_i(m_i))
            % no upper epsilons for rightmost boxes
            inv2 = sprintf('%s %s %s &lt;= x%d &lt;= %s &amp;\n', inv2, ...
                num2str_wrap(partition_i(index)), epsilon_minus, ...
                m_i, num2str_wrap(partition_i(index + 1)));
        else
            % both epsilons for inner boxes
            inv2 = sprintf('%s %s %s &lt;= x%d &lt;= %s %s &amp;\n', inv2, ...
                num2str_wrap(partition_i(index)), epsilon_minus, ...
                m_i, num2str_wrap(partition_i(index + 1)), epsilon_plus);
        end
    end 
    
    % add time partition
    switch (stimulus_mode)
        case 0
            if (G_ADD_TIME)
                inv_string = sprintf('%s 0 &lt;= t &lt; %d &amp;\n %s\n\t\t\t', ...
                                     inv2, time_constraint, inv1);
            else
                inv_string = sprintf('%s%s\n\t\t\t', inv2, inv1);
            end
        
        case {1, 2}
            t_partition = stimulus_wrapper{3};
            current_index = indices(stimulus_dimension);
            if (stimulus_mode == 1)
                lower = num2str_wrap(t_partition(current_index));
                upper = num2str_wrap(t_partition(current_index + 1));
            else
                lower = num2str_wrap(t_partition(current_index + 1));
                upper = num2str_wrap(t_partition(current_index));
            end
            
            if (current_index == 1)
                t_epsilon_str = '';
            else
                t_epsilon_str = epsilon_time;
            end
            
            inv_string = sprintf('%s %s &lt;= t &lt;= %s%s &amp;\n %s\n\t\t\t', ...
                                 inv2, lower, upper, t_epsilon_str, inv1);
        
        otherwise
            assert(false, 'Unknown stimulus mode.');
    end
end


function [ result ] = spaceex_check_property(filename_raw, filepathname_result, is_use_cfg, initial_names, forbidden_names)
% runs SpaceEx on the given model file and writes the result to the given output
% file

    global spaceex_initial_states;
    global spaceex_forbidden_states;
    global time_constraint;
    global G_TRANSITION_EPSILONS;
    global G_ADD_TIME;
    global G_CONFIG_FILE;
    global G_SPACEEX_DEBUG_OUTPUT;
    
    PREFIX = './exports/';
    spaceex_path    = './spaceex/spaceex';
    CFG_NAME        = [PREFIX, filename_raw, '.cfg'];
    spaceex_command = [spaceex_path ' -m ', PREFIX, filename_raw, '.xml'];
    % add screen file
    screen_file     = [PREFIX, filename_raw, '-screen.txt'];
    output_file     = [PREFIX, 'output.txt'];
    output_format      = 'INTV'; % INTV, TXT
    if (is_use_cfg)
        spaceex_command = [spaceex_command, ' -g ', CFG_NAME];
    else
        rel_err            = '1.0E-12';
        abs_err            = '1.0E-15';
        scenario           = 'phaver'; % phaver, supp, stc
        directions         = 'oct'; % box, oct, uni32, ...
        directions2        = ''; % leave empty to ignore
        set_aggregation    = 'none'; % chull, none
        time_horizon       = sprintf('%d', time_constraint + 1);
        sampling_time      = '0.1';
        flowpipe_tolerance = '0.01';
        iter_max           = '1000';
        initial_states     = spaceex_initial_states;
        % add transition epsilons to initial states
        if (G_TRANSITION_EPSILONS > 0)
            initial_states = strcat(initial_states, ' & t_eps==0 & eps==', ...
                num2str(G_TRANSITION_EPSILONS));
        end
        % add discrete locations to initial states
        if (~ isempty(initial_states))
            initial_states = ...
                [initial_states, ' & (loc() == ', initial_names{1}];
        else
            initial_states = ['(loc() == ', initial_names{1}];
        end
        for i = 2 : size(initial_names, 1)
            initial_states = ...
                [initial_states, ' | loc() == ', initial_names{i}];
        end
        initial_states = [initial_states, ')'];
        forbidden_states = spaceex_forbidden_states;
        % add discrete locations to forbidden states
        if (~ isempty(forbidden_states))
            forbidden_states = [forbidden_states, ' & (loc() == '];
        else
            forbidden_states = '(loc() == ';
        end
        or = '';
        for i = 1 : size(forbidden_names, 1)
            if (isempty(forbidden_names{i}))
                continue;
            end
            forbidden_states = [forbidden_states, or, forbidden_names{i}];
            if (isempty(or))
                or = ' | loc() == ';
            end
        end
        forbidden_states = [forbidden_states, ')'];
        % when finding the forbidden states, we can stop
        search_function    = 'reach-verimag';
        % only when an intersection was found, it should be printed.
        output_intersection_only = 'true';
        lha_system         = 'pma';
        if (G_ADD_TIME)
            output_variables   = 't,x1';
        else
            output_variables   = 'x1';
        end

        if (G_CONFIG_FILE)
            % write a config file
            CFG_FILE = fopen(CFG_NAME, 'w');
            spaceex_command = [spaceex_command, ' -g ', CFG_NAME];

            fprintf(CFG_FILE, ['scenario = ', scenario]);
            % write parameters depending on the scenario
            switch (scenario)
                case 'phaver'
                    % PHAVer scenario

                case {'supp', 'stc'}
                    % common settings for LGG and STC scenarios
                    fprintf(CFG_FILE, ['\nrel-err = ', rel_err]);
                    fprintf(CFG_FILE, ['\nabs-err = ', abs_err]);
                    fprintf(CFG_FILE, ['\nset-aggregation = ', set_aggregation]);

                otherwise
                    assert(false, 'Unknown scenario.');
            end

            % write scenario-independent parameters
            fprintf(CFG_FILE, ['\ndirections = ', directions]);
            if (~ isempty(directions2))
                fprintf(CFG_FILE, ['\ndirections = ', directions2]);
            end
            fprintf(CFG_FILE, ['\ntime-horizon = ', time_horizon]);
            fprintf(CFG_FILE, ['\nsampling-time = ', sampling_time]);
            fprintf(CFG_FILE, ['\niter-max = ', iter_max]);
            fprintf(CFG_FILE, ['\nflowpipe-tolerance = ', flowpipe_tolerance]);
            fprintf(CFG_FILE, ['\noutput-format = ', output_format]);
            fprintf(CFG_FILE, ['\noutput-variables = ', output_variables]);
            fprintf(CFG_FILE, ['\noutput-file = ', output_file]);
            fprintf(CFG_FILE, ['\nsystem = ', lha_system]);
            fprintf(CFG_FILE, ['\ninitially = "', initial_states, '"']);
            fprintf(CFG_FILE, ['\nforbidden = "', forbidden_states, '"']);
            fprintf(CFG_FILE, ['\nsearch-function = ', search_function]);
            fprintf(CFG_FILE, ['\noutput-intersection-only = ', ...
                               output_intersection_only]);

            fclose(CFG_FILE);
        else
            % pass the config directly in the command
            spaceex_command = [spaceex_command, ...
                ' --rel-err ''' rel_err ...
                ''' --abs-err ''' abs_err ...
                ''' -c ''' scenario ...
                ''' --directions ''' directions];
            if (~ isempty(directions2))
                % combines both directions settings
                spaceex_command = [spaceex_command, ...
                    ''' --directions ''' directions2];
            end
            spaceex_command = [spaceex_command, ...
                ''' --set-aggregation ''' set_aggregation ...
                ''' --time-horizon ''' time_horizon ...
                ''' --sampling-time ''' sampling_time ...
                ''' --flowpipe-tolerance ''' flowpipe_tolerance ...
                ''' --iter-max ''' iter_max ...
                ''' --system ''' lha_system ...
                ''' -i ''' initial_states ...
                ''' --forbidden ''' forbidden_states ...
                ''' --search-function ''' search_function ...
                ''' --output-intersection-only ''' output_intersection_only ...
                ''' -f ''' output_format ...
                ''' -a ' output_variables ...
                ' -o ' output_file ...
                ];
        end
    end
    
    % add screen file
    % NOTE: We use the file to determine the output of SpaceEx.
    if (G_SPACEEX_DEBUG_OUTPUT)
        % write all debug output
        spaceex_command = [spaceex_command, ' -v D7 > ' screen_file];
    else
        % write only the smallest amount of information
        spaceex_command = [spaceex_command, ' -v l > ' screen_file];
    end
    
    [s, answer] = system(spaceex_command);

    if s~=0 %error
        fprintf(['\n Error when model checking: \n' answer '\n'])
        result = 2;
        return;
    end


    fid=fopen(output_file,'rt');
    if ~feof(fid)
        tline = fgetl(fid);
        % TXT format
        if (strcmp(output_format, 'TXT'))
            tline(1:2)
            if strcmp(tline(1:2), '{}')
                result = 1;
            else
                result = 0;
            end
        % INTV format: file content is 'empty set' for empty intersection
        elseif (strcmp(output_format, 'INTV'))
            if (length(tline) >= 5 && strcmp(tline(1:5), 'empty'))
                [~, line_numbers] = ...
                    system(['grep "Found fixpoint after" ' screen_file ' | wc -l']);
                if (str2double(line_numbers) == 1)
                    % fixed point found, answer is 1 (not excited)
                    result = 1;
                else
                    [~, line_numbers] = ...
                        system(['grep "Performed max. number of iterations" ' screen_file ' | wc -l']);
                    if (str2double(line_numbers) == 1)
                        % no fixed point found, answer is 2 (unknown)
                        result = 2;
                    else
                        fprintf(['\n Error: Strange output in screen file: ' screen_file '\n'])
                        result = 2;
                    end
                end
            else
                % intersection found, answer is 0 (excited)
                result = 0;
            end
        else
            fprintf(['\n Unknown output format: \n' output_format '\n']);
            result = 2;
        end
    else
        fprintf('\n Output file is empty!\n');
        result = 2;
    end
    fclose(fid);

    % printing result
%     output(['result = ' num2str(result) '\n']);
    newfile = fopen(filepathname_result, 'w');
    fprintf(newfile, num2str(result));
    fclose(newfile);
end


function [ result, reachable_locations_out ] = hypy_check_property(filename_raw, filepathname_result, initial_names, reachable_locations_old)
% runs SpaceEx on the given model file and writes the result to the given output
% file

    global spaceex_initial_states;
    global G_FORBIDDEN_LOCS;
    global time_constraint;
    global box_nb_i;
    global G_TRANSITION_EPSILONS;
    global G_ADD_TIME;
    global G_USE_OLD_RESULTS;
    global g_benchmark_reachable_locations_graph_reachability;
    global g_benchmark_reachable_locations_spaceex_reachability_total;
    
    fprintf('HyPy results: ');
    
    PREFIX = ['./exports/', filename_raw];
    
    % write a config file only if initial states are defined
    % (otherwise the model and config files were already existed beforehand)
    if (~ isempty(initial_names))
        % config options
        rel_err            = '1.0E-12';
        abs_err            = '1.0E-15';
        scenario           = 'phaver'; % phaver, supp, stc
        directions         = 'oct'; % box, oct, uni32, ...
        directions2        = ''; % leave empty to ignore
        set_aggregation    = 'none'; % chull, none
        time_horizon       = sprintf('%d', time_constraint + 1);
        sampling_time      = '0.1';
        flowpipe_tolerance = '0.01';
        output_format      = 'INTV'; % INTV, TXT
        iter_max           = '1000';
        lha_system         = 'pma';

        % initial states
        initial_states = spaceex_initial_states;
        % add transition epsilons to initial states
        if (G_TRANSITION_EPSILONS > 0)
            initial_states = strcat(initial_states, ' & t_eps==0 & eps==', ...
                num2str(G_TRANSITION_EPSILONS));
        end
        % add discrete locations to initial states
        if (~ isempty(initial_states))
            initial_states = ...
                [initial_states, ' & (loc() == ', initial_names{1}];
        else
            initial_states = ['(loc() == ', initial_names{1}];
        end
        for i = 2 : size(initial_names, 1)
            initial_states = ...
                [initial_states, ' | loc() == ', initial_names{i}];
        end
        initial_states = [initial_states, ')'];

        if (G_ADD_TIME)
            output_variables = 't,x1';
        else
            output_variables = 'x1';
        end

        % write config file
        cfg_file = fopen([PREFIX, '.cfg'], 'w');
        fprintf(cfg_file, ['scenario = ', scenario]);
        % write parameters depending on the scenario
        switch (scenario)
            case 'phaver'
                % PHAVer scenario

            case {'supp', 'stc'}
                % common settings for LGG and STC scenarios
                fprintf(cfg_file, ['\nrel-err = ', rel_err]);
                fprintf(cfg_file, ['\nabs-err = ', abs_err]);
                fprintf(cfg_file, ['\nset-aggregation = ', set_aggregation]);

            otherwise
                assert(false, 'Unknown scenario.');
        end
        % write scenario-independent parameters
        fprintf(cfg_file, ['\ndirections = ', directions]);
        if (~ isempty(directions2))
            fprintf(cfg_file, ['\ndirections = ', directions2]);
        end
        fprintf(cfg_file, ['\ntime-horizon = ', time_horizon]);
        fprintf(cfg_file, ['\nsampling-time = ', sampling_time]);
        fprintf(cfg_file, ['\niter-max = ', iter_max]);
        fprintf(cfg_file, ['\nflowpipe-tolerance = ', flowpipe_tolerance]);
        fprintf(cfg_file, ['\noutput-format = ', output_format]);
        fprintf(cfg_file, ['\noutput-variables = ', output_variables]);
        fprintf(cfg_file, ['\nsystem = ', lha_system]);
        fprintf(cfg_file, ['\ninitially = "', initial_states, '"']);
        fclose(cfg_file);
    end
    
    % run HyPy (only if no output file exists, i.e., not run beforehand)
    output_file_name = [PREFIX, '.out'];
    if ((~ G_USE_OLD_RESULTS) || (~ exist(output_file_name, 'file')))
        hypy_command = ['python ', './hypy/get_locations/get_locations.py ', ...
            PREFIX, '.xml ', output_file_name];
        [s, answer] = system(hypy_command);

        if (s ~= 0) % error
            fprintf(['\n Error when model checking: \n' answer '\n'])
            result = 2;
            return;
        end
    end
    
    % read reachable locations (overwrite old result)
    out_file = fopen(output_file_name, 'rt');
    is_fixpoint_found = true;
    is_past_first_line = false;
    reachable_locations_new = cell(prod(box_nb_i), 1);
    i_locs = 0;
    factors = cumprod(box_nb_i);
    factors = [1, factors(1 : end-1)];
    while (~ feof(out_file))
        label = fgetl(out_file);
        if (~ is_past_first_line)
            if (strcmp(label, 'fixpoint found'))
                is_fixpoint_found = true;
            elseif (strcmp(label, 'no fixpoint found'))
                is_fixpoint_found = false;
            else
                fprintf(['\n Error when model checking: \n', ...
                    'HyPy could not read fixpoint information of SpaceEx.\n']);
                result = 2;
                return;
            end
            is_past_first_line = true;
            continue;
        end
        i_locs = i_locs + 1;
        idx = getIndexFromLabel(label, factors);
        reachable_locations_new{i_locs} = idx;
    end
    fclose(out_file);
    
    % remove trailing entries and sort reachable locations
    reachable_locations_new(i_locs + 1 : end) = [];
    reachable_locations_new = sortrows(reachable_locations_new);
    
    % check for bad locations
    is_intersection_found = false;
    i_forbidden = 1;
    for i = 1 : length(reachable_locations_new)
        i_reachable = reachable_locations_new{i};
        while ((G_FORBIDDEN_LOCS(i_forbidden) < i_reachable) && ...
                (i_forbidden < length(G_FORBIDDEN_LOCS)))
            i_forbidden = i_forbidden + 1;
        end
        if (G_FORBIDDEN_LOCS(i_forbidden) == i_reachable)
            % a bad location is reachable
            is_intersection_found = true;
            break;
        elseif ((i_forbidden == length(G_FORBIDDEN_LOCS)) && ...
                (G_FORBIDDEN_LOCS(i_forbidden) < i_reachable))
            % checked all bad locations
            break;
        end
    end
    
    if (is_intersection_found)
        % intersection found, answer is 0 (unsafe)
        result = 0;
    elseif (~ is_fixpoint_found)
        % fixpoint was not reached, answer is 2 (unknown)
        result = 2;
    else
        % no intersection found, answer is 1 (safe)
        result = 1;
    end
    
    % writing result file
    res_file = fopen(filepathname_result, 'w');
    fprintf(res_file, num2str(result));
    fclose(res_file);
    
    % printing results
%     output(['result = ' num2str(result) '\n']);
    if (is_fixpoint_found)
        str_fixpoint = '';
        reachable_locations_out = reachable_locations_new;
    else
        str_fixpoint = 'at least ';
        reachable_locations_out = reachable_locations_old;
    end
    if (g_benchmark_reachable_locations_graph_reachability > 0)
        if (result == 0)
            % update benchmark counter only for unsafe cases
            g_benchmark_reachable_locations_spaceex_reachability_total = ...
                g_benchmark_reachable_locations_spaceex_reachability_total + ...
                (g_benchmark_reachable_locations_graph_reachability - length(reachable_locations_new));
        end
        
        output(sprintf('<strong>%s%d/%d locations are reachable (- %d)</strong>\n', ...
            str_fixpoint, length(reachable_locations_new), ...
            g_benchmark_reachable_locations_graph_reachability, ...
            g_benchmark_reachable_locations_graph_reachability - length(reachable_locations_new)));
    else
        output(sprintf('<strong>%s%d/? locations are reachable</strong>\n', ...
            str_fixpoint, length(reachable_locations_new)));
    end
end


function [ index ] = getIndexFromLabel(label, factors)
% converts a location name to an index

    separators = strfind(label, '_');
    index = 1;
    for i = 1 : length(separators)
        if (i == length(separators))
            terminator = length(label);
        else
            terminator = separators(i + 1) - 1;
        end
        substring = label(separators(i) + 1 : terminator);
        num = str2num(substring);
        index = index + factors(i) * (num - 1);
    end
end


function write_transitions(FILE, TRANS_REL, ID_OFFSET_MAP)
% writes the transitions

    global G_TRANSITION_EPSILONS;
    
    % label positioning
    LABEL_SIZE = 64;
    X_HORIZONTAL = -(LABEL_SIZE / 2);
    X_UP = -LABEL_SIZE;
    X_DOWN = 0;
    Y_VERTICAL = -8;
    Y_LEFT = 16;
    Y_RIGHT = -(LABEL_SIZE / 2);
    
    for i = 1 : size(TRANS_REL, 1)
        label_src_id = TRANS_REL(i, 1);
        label_dest_id = TRANS_REL(i, 2);
        
        % adapt location IDs when some of them were removed
        if (ID_OFFSET_MAP ~= -1)
            src = ID_OFFSET_MAP(label_src_id);
            % skip unreachable source locations
            if (src == 0)
                continue;
            end
            dest = ID_OFFSET_MAP(label_dest_id);
            % skip unreachable target locations
            if (dest == 0)
                continue;
            end
        else
            src = label_src_id;
            dest = label_dest_id;
        end
        
        % skip self-loops
        if (src == dest)
            continue;
        end
        
        % write transition header and label
        transition_tag_open(FILE, src, dest);
         % TODO do we need this?
        label2_tags(FILE, sprintf('T%d_%d', label_src_id, label_dest_id));
        
        if (G_TRANSITION_EPSILONS > 0)
            % write epsilon constraint
            guard_tags(FILE, 't_eps &gt;= eps');
            
            % write reset
            assignment_tags(FILE, 't_eps'' == 0');
        end
        
        % write label position % TODO do we need this?
        indices_src  = find_box_indices(src);
        indices_dest = find_box_indices(dest);
        index_diff = indices_dest - indices_src;
        [dim, ~] = find(index_diff ~= 0);
        assert((size(dim, 2) == 1), ...
            'Transitions should only change one dimension');
        if (index_diff(dim) < 0)
            % decrease in dimension
            if (dim > 1) % dimension >1; downward
                X_POS = X_DOWN;
                Y_POS = Y_VERTICAL;
            else % dimension 1; rightward
                X_POS = X_HORIZONTAL;
                Y_POS = Y_RIGHT;
            end
        else
            % increase in dimension
            if (dim > 1) % dimension >1; upward
                X_POS = X_UP;
                Y_POS = Y_VERTICAL;
            else % dimension 1; leftward
                X_POS = X_HORIZONTAL;
                Y_POS = Y_LEFT;
            end
        end
        labelposition_tags(FILE, X_POS, Y_POS, LABEL_SIZE, LABEL_SIZE);
        
        % finalize transition
        transition_tag_close(FILE);
    end
end


function [ indices ] = find_box_indices(ind)
% computes the box number, given the indices array
% TODO do we need this?

    global box_nb_i;
    
    NUM_VAR = size(box_nb_i, 2);
    indices = zeros(NUM_VAR, 1);
    for j = 1 : NUM_VAR-1
        wid        = box_nb_i(j);
        indices(j) = mod (ind-1, wid) + 1;
        ind        = (ind-1) - mod (ind - 1, wid) + 1;
        ind        = ((ind-1) / wid) + 1;
    end
    indices(NUM_VAR) = ind;
end


function [ str ] = num2str_wrap(num)
% converts a number to a string

    % use MATLAB function 'num2str' with precision 10
    str = num2str(num, 10);
end


function output(str)
% output function to quickly deactivate console output

    global G_DEBUG_OUTPUT;
    if (G_DEBUG_OUTPUT)
        fprintf(str);
    end
end

% --- SpaceEx model file writing functions ---

function spaceex_tag_open(file)
    fprintf (file, '<?xml version="1.0" encoding="iso-8859-1"?>\n');
    fprintf (file, '<sspaceex xmlns="http://www-verimag.imag.fr/xml-namespaces/sspaceex" version="0.2" math="SpaceEx">\n');
end

function param_tag(file, name, type, local, d1, d2, dynamics)
    fprintf (file, '\t\t<param name="%s" type="%s" local="%s" d1="%s" d2="%s" dynamics="%s" />\n',name, type, local, d1, d2, dynamics);
end

function param_tag_controlled(file, name, type, local, d1, d2, dynamics, controlled)
    fprintf (file, '\t\t<param name="%s" type="%s" local="%s" d1="%s" d2="%s" dynamics="%s" controlled="%s" />\n',name, type, local, d1, d2, dynamics, controlled);
end

function param_tag_label(file, name, local)
    fprintf (file, '\t\t<param name="%s" type="label" local="%s" />\n', name, local);
end

function flow_tags(file, flow)
    fprintf (file, '\t\t\t<flow>%s</flow>\n', flow);
end

function label2_tags(file, label)
    fprintf (file, '\t\t\t<label>%s</label>\n', label);
end

function labelposition_tags(file, x, y, width, height)
    fprintf (file, '\t\t\t<labelposition x="%d" y="%d" width="%d" height="%d"></labelposition>\n', x, y, width, height);
end

function guard_tags(file, guard)
    fprintf (file, '\t\t\t<guard>%s</guard>\n', guard);
end

function assignment_tags(file, assignment)
    fprintf (file, '\t\t\t<assignment>%s</assignment>\n', assignment);
end

function  invariant_tags(file, invariant)
    fprintf (file, '\t\t\t<invariant>%s</invariant>\n', invariant);
end

function location_tag_open(file, id, name, x, y, width, height)
    fprintf (file, '\t\t<location id="%d" name="%s" x="%d" y="%d"  width="%d" height="%d">\n', id, name, x, y, width, height);
end

function location_tag_close(file)
    fprintf (file, '\t\t</location>\n');
end

function component_tag_open(file, id)
    fprintf (file, '\t<component id="%s">\n', id);
end

function component_tag_close(file)
    fprintf (file, '\t</component>\n');
end
         
function transition_tag_open(file, source, target)
    fprintf (file, '\t\t<transition source="%d" target="%d">\n', source, target);
end

function transition_tag_close(file)
    fprintf (file, '\t\t</transition>\n');
end

function spaceex_tag_close(file)
    fprintf (file, '</sspaceex>\n');
end