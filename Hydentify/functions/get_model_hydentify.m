function  get_model_hydentify()
%get_model - load a PWMA model, TL properties and computation options
%
% Syntax: get_model()
%
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
%
% The information is either given directly by the user, or parsed from a
% file, or loaded from a previously-saved model 
% The user is prompted for choosing between these modes.
% The model is stored in global variables.
%
% -------------------------------------------------------------------------
% Author:  Ezio Bartocci
%          State University of New York at Stony Brook 
% email:   ezio.bartocci@gmail.com
% Website: http://www.eziobartocci.com
%
% Author:  Gregory Batt
%          Boston University, 
%          Brookline, MA, USA
% email:   batt@bu.edu
% Website: http://iasi.bu.edu/~batt/
%
% Author:  Radu Grosu
%          State University of New York at Stony Brook        
% email:   grosu@cs.stonybrook.edu
% Website: http://www.cs.sunysb.edu/~grosu/
%
% Author:  Christian Schilling
%          University of Freiburg, Germany
% email:   schillic at informatik.uni-freiburg.de
% Website: http://swt.informatik.uni-freiburg.de/staff/christian_schilling
%
% October 2014
% -------------------------------------------------------------------------

global analysis_type;

% total time
global time_constraint;

skip = 1;

if (skip ~= 1)
    analysis_type = input('Type of analysis?\n  0: robust property verification,\n  1: search for valid parameter sets.\n 0/1 [0]?: ');
    
    time_constraint = input('Total time constraint? \n');
end

parse_model();
% end of get_model
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% 
% -------------------------------------------------------------------------
function  parse_model()
%parse_model - read input from data file and store information in data structures
%
% Syntax: parse_model() reads input and store information in
% data structures. Information is given as a text file, assumed to be
% present in folder ./models
%
% Called by get_model 
% Calls parse_ramp_expr and read_data (see below)
%

global model_name; 
global variable_nb;
global partition; 
global box_nb_i;
global input_variable_nb;
global is_input_variable; 
global unknown_parameter_nb;
global unknown_parameter; 
global production_rate_parameter;
global degradation_rate_parameter;
global production_ramp_expression;
global degradation_ramp_expression; 
global atomic_proposition; 
global formula;
global init_formula;
global compute_SCC;
global spaceex_initial_states;
global spaceex_forbidden_states;
global stimulus_mode;
global stimulus_wrapper;
global time_constraint;
global G_INITIAL_LOCS;
global G_FORBIDDEN_LOCS;
global G_ADD_TIME;
global G_ANALYSIS_STRATEGY;

% model_name = '...';

% ask for model name if none was specified
if (isempty(model_name))
    model_name = input('Model name: ', 's');
end

fprintf('model name = %s', model_name);

fid = fopen(['models/' model_name '.txt'], 'r');


% some default values
formula              = struct('type', 1, 'text', '');
unknown_parameter_nb = 0;
input_variable_nb    = 0;


% -- reading the number of variables --

% number of (state and input) variables
variable_nb = read_data(fid, 0);


% -- reading the partition --

partition = cell(1, variable_nb);
box_nb_i  = zeros(1, variable_nb);

for i=1:variable_nb
    % row vector with landmarks for gene i (including 0 and max)
    partition{i} = read_data(fid, 0);
    
    % number of intervals for gene i
    box_nb_i(i)  = length(partition{i}) - 1;
end


% -- reading the production and degradation terms and the parameters --
% These are a combination of row vectors for parameter j of gene i
% (eiter [value] or [min max]) and a ramp expression.

production_rate_parameter   = cell(1,variable_nb);
degradation_rate_parameter  = cell(1,variable_nb);
production_ramp_expression  = cell(1,variable_nb);
degradation_ramp_expression = cell(1,variable_nb);

for i=1:variable_nb % for each gene ...
    % -- reading the production terms --
    
    % number of production terms for gene i
    rate_nb = read_data(fid, 0);
    
    production_rate_parameter{i}  = cell(1, rate_nb);
    production_ramp_expression{i} = cell(1, rate_nb);
    
    for j=1:rate_nb % for each rate term ...
        production_rate_parameter{i}{j} = read_data(fid, 0);
        if size(production_rate_parameter{i}{j}, 2) == 2
            unknown_parameter_nb = unknown_parameter_nb + 1;
            unknown_parameter{unknown_parameter_nb} = [1 i j];
        end
        % production ramp expression as a tree
        [production_ramp_expression{i}{j}, ~] = parse_ramp_expression(1, fid);
    end
    
    % -- reading the degradation terms --
    
    % number of degradation terms for gene i
    rate_nb = read_data(fid, 0);
    
    degradation_rate_parameter{i}  = cell(1, rate_nb);
    degradation_ramp_expression{i} = cell(1, rate_nb);
    
    % detection whether we have an 'input variable' or a 'state variable'
    is_input_variable(i) = (rate_nb == 0);
    if (rate_nb == 0)
        input_variable_nb = input_variable_nb + 1;
    end
    
    for j=1:rate_nb % for each rate term ...
        degradation_rate_parameter{i}{j} = read_data(fid, 0);
        if size(degradation_rate_parameter{i}{j}, 2) == 2
            unknown_parameter_nb = unknown_parameter_nb + 1;
            unknown_parameter{unknown_parameter_nb} = [-1 i j];
        end
        % degradation ramp expression as a tree 
        [degradation_ramp_expression{i}{j}, ~] = parse_ramp_expression(1, fid);
    end
end


% -- reading formula related part --

% - reading atomic propositions -

% number of atomic propositions
atomic_proposition_nb = read_data(fid,0);
atomic_proposition = cell(1, atomic_proposition_nb);
for i=1:atomic_proposition_nb
    % atomic proposition prop_i, of type x_i < lambda_i^j or x_i >
    % lambda_i^j, given as {i ''<'' j} or {i ''>'' j}
    atomic_proposition{i} = eval(read_data(fid, 1));
end

% - reading initial conditions -

% either in NuSMV syntax using prop_i for atomic propositions, or All
init_formula = deblank(read_data(fid, 1));

% - reading formula and related options -

% type of temporal logic for formula
% 0: CTL
% 1: LTL
formula.type = read_data(fid, 0);

% LTL formula in NuSMV syntax and using prop_i for atomic propositions
formula.text = read_data(fid, 1);

% computation of strongly connected components?
% 0: No
% 1: Yes
compute_SCC = read_data(fid, 0);


% -- reading new features needed for Hydentify --

% - reading SpaceEx initial states -
spaceex_initial_states = read_data(fid, 1);
assert(~ isempty(spaceex_initial_states), 'The empty string is not allowed.');
% special handlings
if (spaceex_initial_states(1) == '!')
    % special mode for initial and forbidden states: use the same boxes as
    % RoVerGeNe
    %
    % syntax: ! [i1 ... ik]
    % where ij is a proposition index
    % semantics: conjunction prop(i1) & ... & prop(ik)
    %
    % two benefits:
    % 1) avoid modeling errors by only changing the propositions
    % 2) automatically compute the initial locations for graph reachability
    is_compute_spaceex_locations = true;
    bounds = parse_spaceex_bounds_from_props(...
        str2num(spaceex_initial_states(2 : length(spaceex_initial_states))));
    spaceex_initial_states = bounds_to_string(bounds);
elseif (spaceex_initial_states(1) == '?')
    % special mode for initial and forbidden states: use an easy to parse format
    %
    % syntax: ? \n [l1 u1] \n ... \n [lk uk]
    % semantics: lj (uj) is the lower (upper) bound on variable k
    %
    % benefit:
    % - automatically compute the initial locations for graph reachability
    is_compute_spaceex_locations = true;
    bounds = parse_spaceex_bounds_from_intervals(fid);
    spaceex_initial_states = bounds_to_string(bounds);
else
    is_compute_spaceex_locations = false;
end

if (G_ADD_TIME)
    spaceex_initial_states = [spaceex_initial_states, ' & (t == 0)'];
end

% compute initial locations if we want to use Hydentify
if (is_compute_spaceex_locations)
    switch (G_ANALYSIS_STRATEGY)
        case EOptions.ANALYSIS_ROVERGENE
            % do nothing

        case {EOptions.ANALYSIS_SPACEEX_ONLY, EOptions.ANALYSIS_RECURSIVE, ...
              EOptions.ANALYSIS_EAGER, EOptions.ANALYSIS_LAZY, ...
              EOptions.ANALYSIS_PARALLEL, EOptions.ANALYSIS_EAGER_PARALLEL}
            % compute initial locations
            G_INITIAL_LOCS = compute_spaceex_locations(bounds);

        otherwise
            assert(false, 'Unknown strategy.');
    end
end


% - reading SpaceEx forbidden states -
spaceex_forbidden_states = read_data(fid, 1);
% special modes (see initial states)
assert(~ isempty(spaceex_forbidden_states), 'The empty string is not allowed.');
if (spaceex_forbidden_states(1) == '!')
    is_compute_spaceex_locations = true;
    bounds = parse_spaceex_bounds_from_props(...
        str2num(spaceex_forbidden_states(2 : length(spaceex_forbidden_states))));
    spaceex_forbidden_states = bounds_to_string(bounds);
elseif (spaceex_forbidden_states(1) == '?')
    is_compute_spaceex_locations = true;
    bounds = parse_spaceex_bounds_from_intervals(fid);
    spaceex_forbidden_states = bounds_to_string(bounds);
else
    is_compute_spaceex_locations = false;
end

% compute forbidden locations if we want to use Hydentify
if (is_compute_spaceex_locations)
    switch (G_ANALYSIS_STRATEGY)
        case EOptions.ANALYSIS_ROVERGENE
            % do nothing

        case {EOptions.ANALYSIS_SPACEEX_ONLY, EOptions.ANALYSIS_RECURSIVE, ...
              EOptions.ANALYSIS_EAGER, EOptions.ANALYSIS_LAZY, ...
              EOptions.ANALYSIS_PARALLEL, EOptions.ANALYSIS_EAGER_PARALLEL}
            % compute forbidden locations
            G_FORBIDDEN_LOCS = compute_spaceex_locations(bounds);

        otherwise
            assert(false, 'Unknown strategy.');
    end
end


% - reading stimulus related part -

% Christian: This new part reads in the specification of the stimulus (u).
% The first input is the mode of the stimulus. The other inputs depend on
% the mode. Here are the modes:
%
% 0 = no stimulus at all
%
% 1 = multiplicative stimulus already encoded in the input model
%
% 2 = nondeterministic stimulus with partitioning of time and stimulus and
%     additional exact time transitions
%
%  example: t = [t1, t2], u = [u1, u2] means u1 <= u = u2 for 0 <= t <= t1,
%  u2 <= u <= 0 for t1 <= t <= t2, and u = 0 for t > t2
%
%  The stimulus is modelled as an additional input variable. We add the
%  partition by hand "from future to past", because RoVerGeNe wants
%  intervals to be monotonically increasing. This does not matter for the
%  computations, but we add the transitions "from past to future" in the
%  end.

% - handling of different modes -
stimulus_mode = read_data(fid, 0);
switch (stimulus_mode)
    case 0 % no stimulus mode
        % nothing to do
        
    case 1 % multiplicative stimulus encoded in input model
        % read input
        stim_dimension      = read_data(fid, 0); % stimulus dimension
        t_partition         = read_data(fid, 0); % t (time)
        stim_partition      = read_data(fid, 0); % u (stimulus)
        
        stimulus_wrapper    = cell(1, 4);
        stimulus_wrapper{1} = stim_dimension;
        stimulus_wrapper{2} = zeros(1, prod(box_nb_i));
        stimulus_wrapper{3} = t_partition;
        stimulus_wrapper{4} = stim_partition;
        
    case 2 % nondeterministic stimulus
        % read input
        input_dimension = read_data(fid, 0); % input dimension
        t_partition     = read_data(fid, 0); % t (time)
        stim_partition  = read_data(fid, 0); % u (stimulus)

        len = length(stim_partition) + 2;
        % all numbers from 1 to length of stimulus partition
        one_to_n      = zeros(1, len);
        for i=1:len
            one_to_n(i) = i;
        end
        % stimulus intervals inverted
        stim_inverted = zeros(1, len);
        for i=3:len
            stim_inverted(i) = stim_partition(len - i + 1);
        end
        % time intervals inverted
        len = length(t_partition) + 2;
        t_inverted = zeros(1, len);
        t_inverted(1) = time_constraint;
        for i=2:len - 1
            t_inverted(i) = t_partition(len - i);
        end

        % add a new input variable for the stimulus
        variable_nb                       = variable_nb + 1;
        stim_dimension                    = variable_nb;
        partition{stim_dimension}         = stim_inverted;
        box_nb_i(stim_dimension)          = ...
            length(partition{stim_dimension}) - 1;
        input_variable_nb                 = input_variable_nb + 1;
        is_input_variable(stim_dimension) = 1;

        % add a new production rate expression for the input dimension
        new_j = length(production_rate_parameter{input_dimension}) + 1;
        production_rate_parameter{input_dimension}{new_j} = 1;
        ramp    = cell(1, 4);
        ramp{1} = 'r';
        ramp{2} = stim_dimension;
        ramp{3} = one_to_n;
        ramp{4} = stim_inverted;
        production_ramp_expression{input_dimension}{new_j} = ramp;

        % add empty production and degradation functions for stimulus variable
        production_rate_parameter{stim_dimension}   = cell(1, 0);
        production_ramp_expression{stim_dimension}  = cell(1, 0);
        degradation_rate_parameter{stim_dimension}  = cell(1, 0);
        degradation_ramp_expression{stim_dimension} = cell(1, 0);

        % the wrapper contains all stimulus information:
        % 1: stimulus dimension
        % 2: transition data structure
        %    (array working as a map 'src index -> OneOf(dest index, 0)')
        % 3: time partition (inverted)
        % 4: stimulus partition (inverted)
        % 5: dimensions of input ramp % TODO remove after debugging
        stimulus_wrapper = cell(1, 5);
        stimulus_wrapper{1} = stim_dimension;
        stimulus_wrapper{2} = zeros(1, prod(box_nb_i));
        stimulus_wrapper{3} = t_inverted;
        stimulus_wrapper{4} = stim_inverted;
        stimulus_wrapper{5} = [input_dimension, new_j];
        
    otherwise % unknown stimulus mode
        assert(false, 'Unknown stimulus mode.');
end



% close text file
fclose(fid);

% end of parse_model
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% 
% -------------------------------------------------------------------------
function [tree, current_token] = parse_ramp_expression(current_token,fid)
% parses recursively a ramp expression and returns a corresponding tree
% structure
%
% Called by parse_model
%

% Get current token (in 1/r/*/-)
read_token = read_data(fid ,1);
switch read_token
    case '1'
        tree{1}=1; % only one cell, having value 1
        current_token=current_token+1;
    case 'r' 
        tree{1}= 'r'; %corresponds to r+(x_a, theta_a^th1, theta_a^th2) 
        tree{2}= read_data(fid,0); %a, index of variable in ramp expression 
        %Row vector with indices of abscissa of break points for gene i (including 0 and max)
        tree{3}= read_data(fid,0);
        %Row vector with value of ordinate of break points for gene i (between 0 and 1)
        % Christian: changed call argument from 0 to -1 (see function details)
        tree{4}= read_data(fid,-1);
        current_token= current_token+4;
    case '*'
        tree{1} = '*'; %corresponds to ramp_expr1 * ramp_expr2 
        [tree{2}, current_token] = parse_ramp_expression(current_token+1,fid); % ramp_expr1
        [tree{3}, current_token] = parse_ramp_expression(current_token,fid); % ramp_expr2
    case '-'
        tree{1} = '-'; %corresponds to 1 - ramp_expr 
        [tree{2}, current_token] = parse_ramp_expression(current_token+1,fid); % ramp_expr
end
% end of parse_ramp_expr
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% 
% -------------------------------------------------------------------------
function result = read_data(fid, read_string)
% Read data from file or from user.
% Reads a new line in file fid and returns its content.
% Comments are ignored (ie lines starting with %).
% either returns a string or an expression (depending on read_string value)
%
% called by parse_model

current_line = strtrim(fgetl(fid));

%remove comments
while (current_line(1) == '%')
    current_line = strtrim(fgetl(fid));
end

if (read_string == 1)
    % read line as a string
    result = current_line;
else
    % read line as a number array
    result = str2num(current_line);
    
    % Christian: error, assumed being due to the new mode with 2D entries
    %            This is now tried to be fixed.
    if (isempty(result) && (read_string == -1))
        % initialize matrix to 2x2 (we do not know the width in advance)
        result = zeros(2, 2);
        
        inner = strtok(current_line, '{%s}');
        re = regexp(inner, '[\d\.[]]*', 'match');
        i = 1;
        j = 1;
        last = length(re);
        while (i <= last)
            num = re{i};
            if (strcmp(num(1), '['))
                % new format: pair of two numbers
                result(1, j) = str2double(strtok(num, '[%f'));
                i = i + 1;
                num = re(i);
                result(2, j) = str2double(strtok(num, '%f]'));
            else
                % old format: just one number, duplicate it
                result(1, j) = str2double(num);
                result(2, j) = str2double(num);
            end
            i = i + 1;
            j = j + 1;
        end
    end
end

% end of read_data
% -------------------------------------------------------------------------

function [ bounds ] = parse_spaceex_bounds_from_props(props)
% parse special SpaceEx initial/forbidden states
% The states are represented as a list of propositions which are used for
% RoVerGeNe as well.

global atomic_proposition;
global partition;
global variable_nb;

bounds = cell(variable_nb, 2);

% fill with proposition constraints
for i = props
    prop = atomic_proposition{i};
    var = prop{1};
    if (prop{2} == '>')
        op = 1;
    else
        assert(prop{2} == '<', 'Unknown relation.');
        op = 2;
    end
    bound = partition{var}(prop{3});

    bounds{var, op} = bound;
end


function [ bounds ] = parse_spaceex_bounds_from_intervals(fid)
% parse special SpaceEx initial/forbidden states
% The states are represented as a sequence of intervals, each in a new line.

global variable_nb;

bounds = cell(variable_nb, 2);
for var = 1 : variable_nb
    bounds_pair = read_data(fid, 0);
    for i = 1 : 2
        if (bounds_pair(i) == -1)
            bounds{var, i} = [];
        else
            bounds{var, i} = bounds_pair(i);
        end
    end
end


function [ str ] = bounds_to_string(bounds)
% converts a list of bounds to a SpaceEx string

str = '';
amp = '';
for i = 1 : size(bounds, 1)
    if (isempty(num2str(bounds{i, 1})))
        if (isempty(num2str(bounds{i, 2}))) % no bounds, skip
            continue;
        else % no lower bound
            str = [str, amp, '(', 'x', num2str(i), ' <= ', ...
                   num2str(bounds{i, 2}), ')'];
        end
    else
        if (isempty(num2str(bounds{i, 2}))) % no upper bound
            str = [str, amp, '(', num2str(bounds{i, 1}), ' <= x', ...
                num2str(i), ')'];
        elseif (bounds{i, 1} == bounds{i, 2}) % equality
            str = [str, amp, '(x', num2str(i), ' == ', ...
                num2str(bounds{i, 1}), ')'];
        else % distinct lower and upper bounds
            str = [str, amp, '(', num2str(bounds{i, 1}), ' <= x', ...
                num2str(i), ' <= ', num2str(bounds{i, 2}), ')'];
        end
    end

    amp = ' & ';
end


function [ spaceex_locations ] = compute_spaceex_locations(bounds)
% computes the list of SpaceEx locations, given the bounds for each variable

global partition;
global variable_nb;
global box_nb_i;

total_boxes = 1;

% replace empty entries with extreme bounds
for var = 1 : variable_nb
    if (isempty(bounds{var, 1}))
        bounds{var, 1} = partition{var}(1);
    end
    if (isempty(bounds{var, 2}))
        bounds{var, 2} = partition{var}(length(partition{var}));
    end
end

% find bounds for each dimension (in the sense of partition indices)
partition_bounds = zeros(variable_nb, 2);
for var = 1 : variable_nb
    bound = bounds(var, :);
    partition_i = partition{var};
    
    if (isempty(bound{1, 1}))
        lower = 1;
    else
        lower = find_bound(partition_i, bound{1, 1}, -1);
    end
    if (isempty(bound{1, 2}))
        upper = length(partition_i);
    else
        upper = find_bound(partition_i, bound{1, 2}, bound{1, 1});
    end
    
    if (upper == length(partition_i))
        if (lower == upper)
            lower = upper - 1;
        end
    else
        assert((upper - lower) >= 1, 'Invalid box.');
    end
    
    partition_bounds(var, 1) = lower;
    partition_bounds(var, 2) = upper - 1;
    total_boxes = total_boxes * (upper - lower);
end

% construct box indices from partition bounds
indices = ones(1, variable_nb);
for var = 1 : variable_nb
    indices(1, var) = partition_bounds(var, 1);
end
spaceex_locations = zeros(1, total_boxes);
for i = 1 : total_boxes
    if (i > 1)
        % increment indices
        for var = 1 : variable_nb
            if (indices(var) == partition_bounds(var, 2))
                indices(var) = partition_bounds(var, 1);
            else
                indices(var) = indices(var) + 1;
                break;
            end
        end
    end
    % construct corresponding box index
    box_id = indices(1);
    product = 1;
    for j = 2 : variable_nb
        product = product * box_nb_i(j - 1);
        box_id = box_id + (product * (indices(j) - 1));
    end
    
    spaceex_locations(1, i) = box_id;
end


function [ bound ] = find_bound(partition, exact_bound, lower)
% finds the bound in the state space partition    

bound = length(partition);
for i = 2 : bound
    if (exact_bound < partition(i))
        if (lower == exact_bound)
            bound = i;
        else
            bound = i - 1;
        end
        break;
    end
end